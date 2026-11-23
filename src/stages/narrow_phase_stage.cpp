// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "world/world_stages.h"

#include "ccd/shape_ccd.h"
#include "narrow_phase/manifold.h"
#include "narrow_phase/sat.h"
#include "solver/integration.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <variant>

namespace sqp {
namespace {

enum class CastChoice { None, Cast };

using CandidatePose = WorldState::NarrowPhaseData::Pose;
using CandidateCacheEntry = WorldState::NarrowPhaseData::CandidateCacheEntry;

struct CandidateBody {
    const ShapeHandle& shape;
    CandidatePose start;
    const CandidatePose& end;
    const Vec2& linear_velocity;
    const real& angular_velocity;
    const BodyConfiguration& configuration;
    const detail::ShapeQueryView& view;

    CandidateBody(const CandidateCacheEntry& cached, tick_t start_time)
        : shape(cached.shape), end(cached.end), linear_velocity(cached.linear_velocity),
          angular_velocity(cached.angular_velocity), configuration(cached.configuration),
          view(cached.view) {
        start = cached.source;
        if (start_time != cached.island_time)
            predict_body_state(cached.source.position, cached.source.rotation,
                               cached.linear_velocity, cached.angular_velocity, {}, real(0),
                               start_time - cached.island_time, start.position, start.rotation);
    }
};

bool load_candidate(const WorldState& world, BodyHandle body, ShapeHandle shape,
                    CandidateCacheEntry& result) {
    result.valid = false;
    const auto dense_index = world.body_pool.dense_index(body);
    BodyHandle owner;
    if (!world.body_pool.is_valid(body) || !world.shape_pool.is_valid(shape) || !dense_index ||
        !world.body_pool.get_shape_owner(shape, owner) || owner != body ||
        !world.body_pool.get_configuration(body, result.configuration))
        return false;
    result.view = detail::ShapeQueryView(world.shape_pool, shape);
    if (!result.view.is_valid())
        return false;
    const auto versions = world.body_pool.motion_version_column(body.kind());
    if (*dense_index >= versions.size())
        return false;
    const auto read_column = [&](BodyScalarColumn column, real& value) {
        const auto values = world.body_pool.scalar_column(body.kind(), column);
        if (*dense_index >= values.size())
            return false;
        value = values[*dense_index];
        return true;
    };
    result.source.rotation = real(0);
    result.angular_velocity = real(0);
    if (!read_column(BodyScalarColumn::PositionX, result.source.position.x) ||
        !read_column(BodyScalarColumn::PositionY, result.source.position.y) ||
        !read_column(BodyScalarColumn::VelocityX, result.linear_velocity.x) ||
        !read_column(BodyScalarColumn::VelocityY, result.linear_velocity.y))
        return false;
    if (body.kind() == BodyConfigs::Kind::Rigid &&
        (!read_column(BodyScalarColumn::Rotation, result.source.rotation) ||
         !read_column(BodyScalarColumn::AngularVelocity, result.angular_velocity)))
        return false;
    const auto island_id = world.body_pool.island_of(body);
    const Island* island = island_id ? world.island_pool.get(*island_id) : nullptr;
    const bool frozen = island_id &&
                        world.island_pool.frozen_until(*island_id) > world.schedule.time;
    BodyConfigs::RegisterType registration;
    if (!world.body_pool.get_register_type(body, registration) ||
        (registration != BodyConfigs::RegisterType::ImmutableStatic && !island) ||
        (island && island->time > world.interval.target_end && !frozen))
        return false;
    result.body = body;
    result.shape = shape;
    result.motion_version = versions[*dense_index];
    result.island_time = frozen ? world.schedule.time
                                : island ? island->time : world.schedule.time;
    if (registration == BodyConfigs::RegisterType::ImmutableStatic || frozen) {
        result.linear_velocity = {};
        result.angular_velocity = real(0);
    }
    result.end = result.source;
    if (registration != BodyConfigs::RegisterType::ImmutableStatic && !frozen)
        predict_body_state(result.source.position, result.source.rotation, result.linear_velocity,
                           result.angular_velocity, {}, real(0),
                           world.interval.target_end - result.island_time, result.end.position,
                           result.end.rotation);
    result.valid = true;
    return true;
}

CastChoice select_mode(const CandidateBody& body, const CcdSettings& settings) {
    if ((body.shape.kind == ShapeKind::Circle && !settings.enable_circle) ||
        (body.shape.kind == ShapeKind::Polygon && !settings.enable_polygon))
        return CastChoice::None;

    const real motion_threshold = settings.linear_start_threshold_m;
    const real linear_threshold_squared = motion_threshold * motion_threshold;
    const real linear_motion_squared = (body.end.position - body.start.position).length_squared();
    const bool rotational_motion =
        math::abs(body.end.rotation - body.start.rotation) * body.view.sweep_radius() >
        motion_threshold;

    switch (body.configuration.configs.ccd_override) {
    case BodyConfigs::CcdOverride::Disabled:
        return CastChoice::None;
    case BodyConfigs::CcdOverride::Shapecast:
        return CastChoice::Cast;

    case BodyConfigs::CcdOverride::Auto:
        if (rotational_motion || linear_motion_squared > linear_threshold_squared)
            return CastChoice::Cast;
    }
    return CastChoice::None;
}

bool overlaps(const CandidateBody& a, const CandidateBody& b, const CandidatePose& pose_a,
              const CandidatePose& pose_b, ContactManifold* manifold = nullptr) {
    if (!manifold && std::holds_alternative<ConvexPoly>(a.view.geometry()) &&
        std::holds_alternative<ConvexPoly>(b.view.geometry())) {
        PreparedPolygon prepared_a, prepared_b;
        const SatHit hit = sat_depth_query(prepared_a.get(a.view), pose_a.position,
                                           pose_a.rotation, prepared_b.get(b.view),
                                           pose_b.position, pose_b.rotation);
        return hit.hit;
    }
    ContactManifold contact;
    const bool hit =
        build_contacts_from_shape_pair(a.view, pose_a.position, pose_a.rotation, b.view,
                                       pose_b.position, pose_b.rotation, contact);
    if (hit && manifold)
        *manifold = contact;
    return hit;
}

bool axis_aligned_rectangle(const CandidateBody& body, const CandidatePose& start,
                            const CandidatePose& end, Aabb& bounds_start, Aabb& bounds_end) {
    const auto* polygon = std::get_if<ConvexPoly>(&body.view.geometry());
    if (!polygon || polygon->local.size() != 4 ||
        math::abs(end.rotation - start.rotation) > config::EPSILON)
        return false;
    const auto bounds_at = [&](const CandidatePose& pose, Aabb& bounds) {
        Vec2 minimum{std::numeric_limits<real>::infinity(), std::numeric_limits<real>::infinity()};
        Vec2 maximum{-std::numeric_limits<real>::infinity(),
                     -std::numeric_limits<real>::infinity()};
        Vec2 vertices[4];
        for (std::size_t i = 0; i < 4; ++i) {
            vertices[i] = pose.position + math::rotate_vector(polygon->local[i], pose.rotation);
            minimum.x = std::min(minimum.x, vertices[i].x);
            minimum.y = std::min(minimum.y, vertices[i].y);
            maximum.x = std::max(maximum.x, vertices[i].x);
            maximum.y = std::max(maximum.y, vertices[i].y);
        }
        if (!std::isfinite(minimum.x) || !std::isfinite(minimum.y) || !std::isfinite(maximum.x) ||
            !std::isfinite(maximum.y) || maximum.x <= minimum.x || maximum.y <= minimum.y)
            return false;
        bool corners[4]{};
        for (Vec2 vertex : vertices) {
            const int x = math::abs(vertex.x - maximum.x) < real(1e-5) ? 1 : 0;
            const int y = math::abs(vertex.y - maximum.y) < real(1e-5) ? 1 : 0;
            const Vec2 corner{x ? maximum.x : minimum.x, y ? maximum.y : minimum.y};
            if ((vertex - corner).length_squared() > real(1e-10) || corners[x + y * 2])
                return false;
            corners[x + y * 2] = true;
        }
        bounds = {minimum, maximum};
        return true;
    };
    return bounds_at(start, bounds_start) && bounds_at(end, bounds_end);
}

bool circle_center_moves_linearly(const CandidateBody& body, const CandidatePose& start,
                                  const CandidatePose& end) {
    return math::abs(end.rotation - start.rotation) <= config::EPSILON ||
           body.view.local_center.length_squared() <= config::EPSILON * config::EPSILON;
}

Vec2 circle_center_at(const CandidateBody& body, const CandidatePose& pose) {
    return pose.position + math::rotate_vector(body.view.local_center, pose.rotation);
}

} // namespace

tick_t ceil_toi_time(tick_t start, tick_t end, real fraction) {
    const tick_t duration = end - start;
    const long double scaled = static_cast<long double>(std::clamp(fraction, real(0), real(1))) *
                               static_cast<long double>(duration);
    const tick_t offset = scaled >= static_cast<long double>(duration)
                              ? duration
                              : static_cast<tick_t>(std::ceil(scaled));
    return start + offset;
}

bool run_narrow_phase(WorldState& state) {
    auto& output = state.narrow_phase;
    auto& pairs = state.pairs;
    const auto& broad = state.broad_phase;

    if (broad.bodies.size() != broad.shapes.size())
        return stage_error(state, "candidate body and shape columns have inconsistent sizes");
    if (output.candidate_cache.size() < broad.bodies.size())
        output.candidate_cache.resize(broad.bodies.size());
    if (++output.candidate_cache_epoch == 0) {
        for (auto& entry : output.candidate_cache)
            entry.epoch = 0;
        output.candidate_cache_epoch = 1;
    }
    const auto cached_at = [&](u32 index) -> const CandidateCacheEntry& {
        auto& entry = output.candidate_cache[index];
        if (entry.epoch != output.candidate_cache_epoch) {
            entry.epoch = output.candidate_cache_epoch;
            load_candidate(state, broad.bodies[index], broad.shapes[index], entry);
        }
        return entry;
    };

    for (const BroadPhase::Pair candidate : broad.candidates) {
        if (candidate.a >= broad.bodies.size() || candidate.b >= broad.bodies.size() ||
            candidate.a == candidate.b)
            continue;

        const BodyHandle body_a = broad.bodies[candidate.a];
        const BodyHandle body_b = broad.bodies[candidate.b];
        if (body_a == body_b)
            continue;

        const auto& cached_a = cached_at(candidate.a);
        const auto& cached_b = cached_at(candidate.b);
        if (!cached_a.valid || !cached_b.valid)
            continue;
        const auto id_a = state.body_pool.island_of(body_a);
        const auto id_b = state.body_pool.island_of(body_b);
        const bool frozen_a = id_a && state.island_pool.frozen_until(*id_a) > state.schedule.time;
        const bool frozen_b = id_b && state.island_pool.frozen_until(*id_b) > state.schedule.time;
        if ((frozen_a && frozen_b) ||
            (frozen_a && cached_b.configuration.configs.register_type ==
                             BodyConfigs::RegisterType::ImmutableStatic) ||
            (frozen_b && cached_a.configuration.configs.register_type ==
                             BodyConfigs::RegisterType::ImmutableStatic))
            continue;

        const tick_t start_time =
            std::max({cached_a.island_time, cached_b.island_time, state.schedule.time});
        CandidateBody a(cached_a, start_time);
        CandidateBody b(cached_b, start_time);
        const Vec2 origin = a.start.position;
        const CandidatePose local_start_a{{}, a.start.rotation};
        const CandidatePose local_end_a{a.end.position - origin, a.end.rotation};
        const CandidatePose local_start_b{b.start.position - origin, b.start.rotation};
        const CandidatePose local_end_b{b.end.position - origin, b.end.rotation};

        const CastChoice mode_a = select_mode(a, state.settings.ccd);
        const CastChoice mode_b = select_mode(b, state.settings.ccd);

        const Pair pair{
            body_a, body_b, a.shape, b.shape, cached_a.motion_version, cached_b.motion_version};
        if (mode_a == CastChoice::None && mode_b == CastChoice::None) {
            ContactManifold manifold;
            if (overlaps(a, b, a.end, b.end, &manifold)) {
                pairs.pairs->push_back(pair);
                if (!pairs.manifolds->append(manifold)) {
                    pairs.pairs->pop_back();
                    return stage_error(state, "run_narrow_phase: check failed (!pairs.manifolds->append(manifold))");
                }
            }
            continue;
        }

        real fraction;
        Vec2 normal;
        const bool box_a = body_a.kind() == BodyConfigs::Kind::AabbBox;
        const bool box_b = body_b.kind() == BodyConfigs::Kind::AabbBox;
        const bool circle_a = body_a.kind() == BodyConfigs::Kind::LiteCircle;
        const bool circle_b = body_b.kind() == BodyConfigs::Kind::LiteCircle;
        const bool circle_shape_a = std::holds_alternative<Circle>(a.view.geometry());
        const bool circle_shape_b = std::holds_alternative<Circle>(b.view.geometry());
        Aabb rectangle_a0, rectangle_a1, rectangle_b0, rectangle_b1;
        const bool rectangle_a =
            !box_a && !circle_a && box_b &&
            axis_aligned_rectangle(a, local_start_a, local_end_a, rectangle_a0, rectangle_a1);
        const bool rectangle_b =
            !box_b && !circle_b && box_a &&
            axis_aligned_rectangle(b, local_start_b, local_end_b, rectangle_b0, rectangle_b1);
        const auto box_bounds = [](Vec2 position, Vec2 half_extents) {
            return Aabb{position - half_extents, position + half_extents};
        };
        bool hit = false;
        if ((box_a && box_b) || (box_a && circle_b) || (circle_a && box_b) || rectangle_a ||
            rectangle_b) {
            BodyProperties properties_a, properties_b;
            if (!state.body_pool.get_properties(body_a, properties_a) ||
                !state.body_pool.get_properties(body_b, properties_b))
                return stage_error(state, "run_narrow_phase: check failed (!state.body_pool.get_properties(body_a, properties_a) || !state.body_pool.get_properties(body_b, properties_b))");
            if (box_a && box_b) {
                hit = sweep_aabb_vs_aabb_toi(
                    box_bounds(local_start_a.position, properties_a.half_extents),
                    box_bounds(local_end_a.position, properties_a.half_extents),
                    box_bounds(local_start_b.position, properties_b.half_extents),
                    box_bounds(local_end_b.position, properties_b.half_extents), true, real(0),
                    fraction, normal);
            } else if (rectangle_a) {
                hit = sweep_aabb_vs_aabb_toi(
                    rectangle_a0, rectangle_a1,
                    box_bounds(local_start_b.position, properties_b.half_extents),
                    box_bounds(local_end_b.position, properties_b.half_extents), true, real(0),
                    fraction, normal);
            } else if (rectangle_b) {
                hit = sweep_aabb_vs_aabb_toi(
                    box_bounds(local_start_a.position, properties_a.half_extents),
                    box_bounds(local_end_a.position, properties_a.half_extents), rectangle_b0,
                    rectangle_b1, true, real(0), fraction, normal);
            } else if (circle_a) {
                hit = sweep_circle_vs_aabb_toi(
                    local_start_a.position, local_end_a.position, properties_a.radius,
                    box_bounds(local_start_b.position, properties_b.half_extents),
                    box_bounds(local_end_b.position, properties_b.half_extents), true, real(0),
                    fraction, normal);
            } else {
                hit = sweep_circle_vs_aabb_toi(
                    local_start_b.position, local_end_b.position, properties_b.radius,
                    box_bounds(local_start_a.position, properties_a.half_extents),
                    box_bounds(local_end_a.position, properties_a.half_extents), true, real(0),
                    fraction, normal);
                normal = -normal;
            }
        } else if (circle_shape_a && circle_shape_b &&
                   circle_center_moves_linearly(a, local_start_a, local_end_a) &&
                   circle_center_moves_linearly(b, local_start_b, local_end_b)) {
            hit = sweep_circle_vs_circle_toi(
                circle_center_at(a, local_start_a), circle_center_at(a, local_end_a),
                std::get<Circle>(a.view.geometry()).radius, circle_center_at(b, local_start_b),
                circle_center_at(b, local_end_b), std::get<Circle>(b.view.geometry()).radius, true,
                real(0), fraction, normal);
        } else if (circle_a && std::holds_alternative<ConvexPoly>(b.view.geometry()) &&
                   math::abs(local_end_b.rotation - local_start_b.rotation) <= config::EPSILON) {
            Vec2 point;
            hit = sweep_circle_vs_poly_toi(
                circle_center_at(a, local_start_a), circle_center_at(a, local_end_a),
                std::get<Circle>(a.view.geometry()).radius, std::get<ConvexPoly>(b.view.geometry()),
                local_start_b.position, local_start_b.rotation, local_end_b.rotation,
                local_end_b.position - local_start_b.position, true, real(0), fraction, normal,
                point);
        } else if (circle_b && std::holds_alternative<ConvexPoly>(a.view.geometry()) &&
                   math::abs(local_end_a.rotation - local_start_a.rotation) <= config::EPSILON) {
            Vec2 point;
            hit = sweep_circle_vs_poly_toi(
                circle_center_at(b, local_start_b), circle_center_at(b, local_end_b),
                std::get<Circle>(b.view.geometry()).radius, std::get<ConvexPoly>(a.view.geometry()),
                local_start_a.position, local_start_a.rotation, local_end_a.rotation,
                local_end_a.position - local_start_a.position, true, real(0), fraction, normal,
                point);
            normal = -normal;
        } else {
            hit = sweep_shape_pair_toi(
                a.view, local_start_a.position, local_start_a.rotation, local_end_a.position,
                local_end_a.rotation, b.view, local_start_b.position, local_start_b.rotation,
                local_end_b.position, local_end_b.rotation, true, real(0), fraction, normal);
        }
        if (!hit)
            continue;

        if (!std::isfinite(fraction) || !std::isfinite(normal.x) || !std::isfinite(normal.y))
            continue;

        if (fraction == real(0) && overlaps(a, b, local_start_a, local_start_b)) {
            const auto kinematics = ccd_contact_kinematics(
                a.view, local_start_a.position, local_start_a.rotation, a.linear_velocity,
                a.angular_velocity, b.view, local_start_b.position, local_start_b.rotation,
                b.linear_velocity, b.angular_velocity, normal);

            if (kinematics.normal_speed > real(0)) {
                const bool reentry = sweep_shape_pair_reentry_toi(
                    a.view, local_start_a.position, local_start_a.rotation, local_end_a.position,
                    local_end_a.rotation, b.view, local_start_b.position, local_start_b.rotation,
                    local_end_b.position, local_end_b.rotation, true, real(0), fraction, normal);

                if (!reentry)
                    continue;
            }
        }

        if (!std::isfinite(fraction) || !std::isfinite(normal.x) || !std::isfinite(normal.y))
            continue;

        const auto impact_pose = [fraction](const CandidatePose& start, const CandidatePose& end) {
            return CandidatePose{start.position + (end.position - start.position) * fraction,
                                 start.rotation + (end.rotation - start.rotation) * fraction};
        };
        const auto impact_a = impact_pose(local_start_a, local_end_a);
        const auto impact_b = impact_pose(local_start_b, local_end_b);

        output.toi_pairs.push_back(pair);
        output.toi_times.push_back(ceil_toi_time(start_time, state.interval.target_end, fraction));
        output.toi_normals.push_back(normal);
        output.toi_impacts.push_back({origin, impact_a.position, impact_b.position,
                                      impact_a.rotation, impact_b.rotation, true});
    }

    state.broad_phase.candidates.clear();

    return true;
}

} // namespace sqp
