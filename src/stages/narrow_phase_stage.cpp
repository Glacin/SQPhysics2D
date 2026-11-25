// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "world/world_stages.h"

#include "ccd/shape_ccd.h"
#include "narrow_phase/manifold.h"
#include "narrow_phase/sat.h"
#include "solver/integration.h"

#include <algorithm>
#include <cmath>
#include <cstring>
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

    CandidateBody(const CandidateCacheEntry& cached, const CandidatePose& start_pose)
        : shape(cached.shape), start(start_pose), end(cached.end),
          linear_velocity(cached.linear_velocity), angular_velocity(cached.angular_velocity),
          configuration(cached.configuration), view(cached.view) {}
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
    const bool frozen =
        island_id && world.island_pool.frozen_until(*island_id) > world.schedule.time;
    const auto registration = result.configuration.configs.register_type;
    if ((registration != BodyConfigs::RegisterType::ImmutableStatic && !island) ||
        (island && island->time > world.interval.target_end && !frozen))
        return false;
    result.body = body;
    result.shape = shape;
    result.motion_version = versions[*dense_index];
    result.island_time = frozen ? world.schedule.time : island ? island->time : world.schedule.time;
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
    result.frozen = frozen;
    result.predicted_start = result.source;
    result.predicted_start_time = result.island_time;
    result.sweep_radius = result.view.sweep_radius();
    result.dimensions_loaded = false;
    result.valid = true;
    return true;
}

CastChoice select_mode(const CandidateBody& body, const CcdSettings& settings,
                       real linear_threshold_squared, real sweep_radius) {
    if ((body.shape.kind == ShapeKind::Circle && !settings.enable_circle) ||
        (body.shape.kind == ShapeKind::Polygon && !settings.enable_polygon))
        return CastChoice::None;

    switch (body.configuration.configs.ccd_override) {
    case BodyConfigs::CcdOverride::Disabled:
        return CastChoice::None;
    case BodyConfigs::CcdOverride::Shapecast:
        return CastChoice::Cast;

    case BodyConfigs::CcdOverride::Auto:
        return math::abs(body.end.rotation - body.start.rotation) * sweep_radius >
                           settings.linear_start_threshold_m ||
                       (body.end.position - body.start.position).length_squared() >
                           linear_threshold_squared
                   ? CastChoice::Cast
                   : CastChoice::None;
    }
    return CastChoice::None;
}

bool overlaps(const CandidateBody& a, const CandidateBody& b, const CandidatePose& pose_a,
              const CandidatePose& pose_b, ContactManifold* manifold = nullptr) {
    if (!manifold && std::holds_alternative<ConvexPoly>(a.view.geometry()) &&
        std::holds_alternative<ConvexPoly>(b.view.geometry())) {
        PreparedPolygon prepared_a, prepared_b;
        const SatHit hit =
            sat_depth_query(prepared_a.get(a.view), pose_a.position, pose_a.rotation,
                            prepared_b.get(b.view), pose_b.position, pose_b.rotation);
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

// Prove separation on a fixed axis for the entire affine translation/rotation sweep.
bool polygon_sweep_is_separated(const CandidateBody& fixed, const CandidateBody& moving,
                                Vec2 outward_normal) {
    const auto* fixed_polygon = std::get_if<ConvexPoly>(&fixed.view.geometry());
    const auto* moving_polygon = std::get_if<ConvexPoly>(&moving.view.geometry());
    if (!fixed_polygon || !moving_polygon || fixed_polygon->local.size() < 3 ||
        moving_polygon->local.size() < 3)
        return false;

    const Vec2d normal(outward_normal);
    const double normal_length = normal.length();
    if (!std::isfinite(normal_length) || normal_length <= config::EPSILON)
        return false;
    const Vec2d axis = normal / normal_length;
    const Vec2d fixed_center(fixed_polygon->center_offset());
    const double fixed_cos = std::cos(static_cast<double>(fixed.start.rotation));
    const double fixed_sin = std::sin(static_cast<double>(fixed.start.rotation));
    double support = -std::numeric_limits<double>::infinity();
    double fixed_scale = 0;
    for (Vec2 vertex : fixed_polygon->local) {
        const Vec2d v(vertex);
        const auto point = math::rotate_vector(v + fixed_center, fixed_cos, fixed_sin);
        const double projection = math::dot(axis, point);
        if (!std::isfinite(projection))
            return false;
        support = std::max(support, projection);
        fixed_scale = std::max(fixed_scale, v.length());
    }

    const Vec2d start = Vec2d(moving.start.position) - Vec2d(fixed.start.position);
    const Vec2d end = Vec2d(moving.end.position) - Vec2d(fixed.start.position);
    const Vec2d translation = end - start;
    const double rotation = moving.start.rotation;
    const double angular_delta = static_cast<double>(moving.end.rotation) - rotation;
    const double moving_cos = std::cos(rotation);
    const double moving_sin = std::sin(rotation);
    const Vec2d moving_center(moving_polygon->center_offset());
    const double scale = 1.0 + fixed_scale + fixed_center.length() + start.length() + end.length() +
                         moving_center.length();
    constexpr double rounding = 64.0 * std::numeric_limits<double>::epsilon();

    for (Vec2 vertex : moving_polygon->local) {
        const Vec2d v(vertex);
        const Vec2d local_point = v + moving_center;
        const Vec2d offset = math::rotate_vector(local_point, moving_cos, moving_sin);
        const double gap = math::dot(axis, start + offset) - support;
        const Vec2d tangent{-offset.y, offset.x};
        const double slope = math::dot(axis, translation + angular_delta * tangent);
        const double curvature_bound = local_point.length() * angular_delta * angular_delta;
        const double end_lower_bound = gap + slope - 0.5 * curvature_bound;
        const double margin = rounding * (scale + v.length() + std::abs(slope) + curvature_bound);
        // The concave quadratic lower bound is smallest at one of the two endpoints.
        if (!std::isfinite(gap) || !std::isfinite(end_lower_bound) ||
            !(gap > margin && end_lower_bound > margin))
            return false;
    }
    return true;
}

bool circle_center_moves_linearly(const CandidateBody& body, const CandidatePose& start,
                                  const CandidatePose& end) {
    return math::abs(end.rotation - start.rotation) <= config::EPSILON ||
           body.view.local_center.length_squared() <= config::EPSILON * config::EPSILON;
}

Vec2 circle_center_at(const CandidateBody& body, const CandidatePose& pose) {
    return pose.position + math::rotate_vector(body.view.local_center, pose.rotation);
}

// Only ordinary boundaries can leave a newly closing circle already inside a wall.
void defer_wall_contacts(WorldState& state, std::size_t first) {
    const tick_t now = state.schedule.time;
    if (state.schedule.last_ordinary_end != now || now >= state.interval.target_end)
        return;

    auto& output = state.narrow_phase;
    for (std::size_t i = first; i < output.toi_times.size(); ++i) {
        // A positive sweep fraction is always rounded to a later time.
        if (output.toi_times[i] != now)
            continue;
        const Pair& pair = output.toi_pairs[i];
        const bool circle_a = pair.body_a.kind() == BodyConfigs::Kind::LiteCircle;
        const bool circle_b = pair.body_b.kind() == BodyConfigs::Kind::LiteCircle;
        if (circle_a == circle_b)
            continue;
        const BodyHandle circle = circle_a ? pair.body_a : pair.body_b;
        const BodyHandle wall = circle_a ? pair.body_b : pair.body_a;
        const ShapeHandle wall_shape = circle_a ? pair.shape_b : pair.shape_a;
        BodyConfigs::RegisterType registration;
        if (wall_shape.kind != ShapeKind::Polygon ||
            !state.body_pool.get_register_type(wall, registration) ||
            registration != BodyConfigs::RegisterType::ImmutableStatic)
            continue;
        const auto id = state.body_pool.island_of(circle);
        const Island* island = id ? state.island_pool.get(*id) : nullptr;
        BodyState motion;
        if (!island || island->time != now || state.island_pool.frozen_until(*id) > now ||
            !state.body_pool.get_state(circle, motion))
            continue;

        const detail::ShapeQueryView a(state.shape_pool, pair.shape_a);
        const detail::ShapeQueryView b(state.shape_pool, pair.shape_b);
        const auto& impact = output.toi_impacts[i];
        ContactManifold manifold;
        if (!build_contacts_from_shape_pair(a, impact.position_a, impact.rotation_a, b,
                                            impact.position_b, impact.rotation_b, manifold))
            continue;
        const Vec2 relative_velocity = circle_a ? -motion.linear_velocity : motion.linear_velocity;
        if (math::dot(relative_velocity, manifold.normal) < -config::EPSILON &&
            std::any_of(manifold.points, manifold.points + manifold.count,
                        [](const ContactPoint& point) { return point.depth > config::EPSILON; })) {
            output.toi_times[i] = now + 1;
            output.toi_normals[i] = -manifold.normal;
        }
    }
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

template <bool cache_results> static bool run_narrow_phase_candidates(WorldState& state) {
    auto& output = state.narrow_phase;
    auto& pairs = state.pairs;
    auto& broad = state.broad_phase;
    const bool freeze_retry = broad.freeze_retry;
    auto& previous_results = output.pair_results;
    auto& next_results = output.next_pair_results;
    if (!broad.reprepare_after_toi || output.pair_result_end != state.interval.target_end)
        previous_results.clear();
    next_results.clear();
    std::size_t previous_result = 0;
    const bool reuse = freeze_retry && output.prepared_time == state.schedule.time &&
                       output.prepared_end == state.interval.target_end &&
                       output.ordinary_contact_keys.size() == pairs.pairs->size();

    if (freeze_retry) {
        // No TOI is selected during a retry, so the delayed columns can hold the prior rows.
        if (reuse) {
            pairs.pairs.swap(pairs.delayed_pairs);
            pairs.manifolds.swap(pairs.delayed_manifolds);
            output.ordinary_contact_keys.swap(output.previous_contact_keys);
        }
        pairs.pairs->clear();
        pairs.manifolds->clear();
        pairs.toi_local_contacts.clear();
        pairs.original_manifolds.clear();
        pairs.original_manifolds_valid = false;
    }
    output.ordinary_contact_keys.clear();

    if (broad.bodies.size() != broad.shapes.size())
        return stage_error(state, "candidate body and shape columns have inconsistent sizes");
    if (output.candidate_cache.size() < broad.bodies.size())
        output.candidate_cache.resize(broad.bodies.size());
    if (++output.candidate_cache_epoch == 0) {
        for (auto& entry : output.candidate_cache)
            entry.epoch = 0;
        output.candidate_cache_epoch = 1;
    }
    const auto cached_at = [&](u32 index) -> CandidateCacheEntry& {
        auto& entry = output.candidate_cache[index];
        if (entry.epoch != output.candidate_cache_epoch) {
            entry.epoch = output.candidate_cache_epoch;
            load_candidate(state, broad.bodies[index], broad.shapes[index], entry);
        }
        return entry;
    };

    const auto& ccd = state.settings.ccd;
    const real linear_threshold_squared =
        ccd.linear_start_threshold_m * ccd.linear_start_threshold_m;
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
    u64 dimensions_requested = 0, dimensions_loaded = 0;
    u64 starts_predicted = 0, predictions_reused = 0;
#endif
    const auto start_at = [&](CandidateCacheEntry& entry, tick_t time) -> const CandidatePose& {
        if (time == entry.island_time)
            return entry.source;
        if (time != entry.predicted_start_time) {
            predict_body_state(entry.source.position, entry.source.rotation, entry.linear_velocity,
                               entry.angular_velocity, {}, real(0), time - entry.island_time,
                               entry.predicted_start.position, entry.predicted_start.rotation);
            entry.predicted_start_time = time;
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
            ++starts_predicted;
#endif
        } else {
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
            ++predictions_reused;
#endif
        }
        return entry.predicted_start;
    };
    const auto load_dimensions = [&](CandidateCacheEntry& entry) {
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
        ++dimensions_requested;
#endif
        if (entry.dimensions_loaded)
            return true;
        BodyProperties properties;
        if (!state.body_pool.get_properties(entry.body, properties))
            return false;
        entry.radius = properties.radius;
        entry.half_extents = properties.half_extents;
        entry.dimensions_loaded = true;
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
        ++dimensions_loaded;
#endif
        return true;
    };

    const auto append_ordinary = [&](BroadPhase::Pair key, const Pair& pair,
                                     const ContactManifold& manifold) {
        if (!pairs.manifolds->append(manifold))
            return false;
        pairs.pairs->push_back(pair);
        output.ordinary_contact_keys.push_back(key);
        return true;
    };
    const auto was_updated = [&](u32 index) {
        return index < broad.selected.size() && broad.selected[index] == broad.selection_epoch;
    };
    const auto append_toi = [&](const Pair& pair,
                                const WorldState::NarrowPhaseData::SweepHit& hit) {
        output.toi_pairs.push_back(pair);
        output.toi_times.push_back(hit.time);
        output.toi_normals.push_back(hit.normal);
        output.toi_impacts.push_back(hit.impact);
    };
    const std::size_t first_toi = output.toi_times.size();
    std::size_t previous_row = 0;
    for (const BroadPhase::Pair candidate : broad.candidates) {
        if (reuse && !was_updated(candidate.a) && !was_updated(candidate.b)) {
            const auto& keys = output.previous_contact_keys;
            while (previous_row < keys.size() && keys[previous_row] < candidate)
                ++previous_row;
            if (previous_row < keys.size() && keys[previous_row] == candidate &&
                pairs.delayed_manifolds->counts[previous_row] != 0 &&
                !append_ordinary(candidate, (*pairs.delayed_pairs)[previous_row],
                                 pairs.delayed_manifolds->at(previous_row)))
                return stage_error(state, "run_narrow_phase: failed to retain ordinary contact");
            continue;
        }
        if (candidate.a >= broad.bodies.size() || candidate.b >= broad.bodies.size() ||
            candidate.a == candidate.b)
            continue;

        const BodyHandle body_a = broad.bodies[candidate.a];
        const BodyHandle body_b = broad.bodies[candidate.b];
        if (body_a == body_b)
            continue;

        auto& cached_a = cached_at(candidate.a);
        auto& cached_b = cached_at(candidate.b);
        if (!cached_a.valid || !cached_b.valid)
            continue;
        const bool frozen_a = cached_a.frozen;
        const bool frozen_b = cached_b.frozen;
        if ((frozen_a && frozen_b) ||
            (frozen_a && cached_b.configuration.configs.register_type ==
                             BodyConfigs::RegisterType::ImmutableStatic) ||
            (frozen_b && cached_a.configuration.configs.register_type ==
                             BodyConfigs::RegisterType::ImmutableStatic))
            continue;

        const tick_t start_time =
            std::max({cached_a.island_time, cached_b.island_time, state.schedule.time});
        CandidateBody a(cached_a, start_at(cached_a, start_time));
        CandidateBody b(cached_b, start_at(cached_b, start_time));
        
        const Vec2 origin = a.start.position;
        const CandidatePose local_start_a{{}, a.start.rotation};
        const CandidatePose local_end_a{a.end.position - origin, a.end.rotation};
        const CandidatePose local_start_b{b.start.position - origin, b.start.rotation};
        const CandidatePose local_end_b{b.end.position - origin, b.end.rotation};

        const CastChoice mode_a =
            select_mode(a, ccd, linear_threshold_squared, cached_a.sweep_radius);
        const CastChoice mode_b =
            select_mode(b, ccd, linear_threshold_squared, cached_b.sweep_radius);

        const Pair pair{
            body_a, body_b, a.shape, b.shape, cached_a.motion_version, cached_b.motion_version};
        const bool use_cast = mode_a == CastChoice::Cast || mode_b == CastChoice::Cast;
        WorldState::NarrowPhaseData::PairResult* row = nullptr;
        constexpr auto immutable = BodyConfigs::RegisterType::ImmutableStatic;
        // Analytic LiteCircle and AABB casts are cheaper than recording their inputs.
        if (cache_results &&
            (!use_cast ||
             (body_a.kind() == BodyConfigs::Kind::Rigid && a.shape.kind == ShapeKind::Polygon &&
              a.configuration.configs.register_type != immutable) ||
             (body_b.kind() == BodyConfigs::Kind::Rigid && b.shape.kind == ShapeKind::Polygon &&
              b.configuration.configs.register_type != immutable))) {
            auto& entry = next_results.emplace_back();
            entry.candidate = candidate;
            entry.pair = pair;
            entry.modes = u8(mode_a == CastChoice::Cast) | (u8(mode_b == CastChoice::Cast) << 1);
            entry.frozen = u8(frozen_a) | (u8(frozen_b) << 1);
            if (use_cast)
                entry.start_time = start_time;
            if (use_cast) {
                entry.motion = {a.end.position.x,    a.end.position.y,    a.end.rotation,
                                b.end.position.x,    b.end.position.y,    b.end.rotation,
                                a.start.position.x,  a.start.position.y,  a.start.rotation,
                                b.start.position.x,  b.start.position.y,  b.start.rotation,
                                a.linear_velocity.x, a.linear_velocity.y, a.angular_velocity,
                                b.linear_velocity.x, b.linear_velocity.y, b.angular_velocity};
            } else {
                entry.motion = {a.end.position.x, a.end.position.y, a.end.rotation,
                                b.end.position.x, b.end.position.y, b.end.rotation};
            }

            while (previous_result < previous_results.size() &&
                   previous_results[previous_result].candidate < candidate)
                ++previous_result;
            if (previous_result < previous_results.size()) {
                const auto& prior = previous_results[previous_result];
                if (prior.candidate == candidate && prior.pair == pair &&
                    prior.start_time == entry.start_time && prior.modes == entry.modes &&
                    prior.frozen == entry.frozen &&
                    std::memcmp(prior.motion.data(), entry.motion.data(), sizeof(entry.motion)) ==
                        0) {
                    entry.result = prior.result;
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
                    ++state.stats.narrow_reused_results;
#endif
                    if (const auto* hit =
                            std::get_if<WorldState::NarrowPhaseData::SweepHit>(&entry.result))
                        append_toi(pair, *hit);
                    if (const auto* manifold = std::get_if<ContactManifold>(&entry.result);
                        manifold && !append_ordinary(candidate, pair, *manifold))
                        return stage_error(state,
                                           "run_narrow_phase: failed to reuse ordinary contact");
                    continue;
                }
            }
            row = &entry;
        }
        if (!use_cast) {
            ContactManifold manifold;
            if (overlaps(a, b, a.end, b.end, &manifold)) {
                if (row)
                    row->result = manifold;
                if (!append_ordinary(candidate, pair, manifold))
                    return stage_error(state,
                                       "run_narrow_phase: failed to append ordinary contact");
            }
            continue;
        }

#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
        state.stats.freeze_retry_sweeps += freeze_retry;
        ++state.stats.narrow_sweeps;
#endif
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
            if (!load_dimensions(cached_a) || !load_dimensions(cached_b))
                return stage_error(state,
                                   "run_narrow_phase: failed to load analytic sweep dimensions");
            if (box_a && box_b) {
                hit = sweep_aabb_vs_aabb_toi(
                    box_bounds(local_start_a.position, cached_a.half_extents),
                    box_bounds(local_end_a.position, cached_a.half_extents),
                    box_bounds(local_start_b.position, cached_b.half_extents),
                    box_bounds(local_end_b.position, cached_b.half_extents), true, real(0),
                    fraction, normal);
            } else if (rectangle_a) {
                hit = sweep_aabb_vs_aabb_toi(
                    rectangle_a0, rectangle_a1,
                    box_bounds(local_start_b.position, cached_b.half_extents),
                    box_bounds(local_end_b.position, cached_b.half_extents), true, real(0),
                    fraction, normal);
            } else if (rectangle_b) {
                hit = sweep_aabb_vs_aabb_toi(
                    box_bounds(local_start_a.position, cached_a.half_extents),
                    box_bounds(local_end_a.position, cached_a.half_extents), rectangle_b0,
                    rectangle_b1, true, real(0), fraction, normal);
            } else if (circle_a) {
                hit = sweep_circle_vs_aabb_toi(
                    local_start_a.position, local_end_a.position, cached_a.radius,
                    box_bounds(local_start_b.position, cached_b.half_extents),
                    box_bounds(local_end_b.position, cached_b.half_extents), true, real(0),
                    fraction, normal);
            } else {
                hit = sweep_circle_vs_aabb_toi(
                    local_start_b.position, local_end_b.position, cached_b.radius,
                    box_bounds(local_start_a.position, cached_a.half_extents),
                    box_bounds(local_end_a.position, cached_a.half_extents), true, real(0),
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
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
                state.stats.freeze_retry_sweeps += freeze_retry;
                ++state.stats.narrow_sweeps;
#endif
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

        const tick_t toi_time = ceil_toi_time(start_time, state.interval.target_end, fraction);
        if (fraction > real(0) && toi_time - start_time == 1) {
            constexpr auto fixed = BodyConfigs::RegisterType::ImmutableStatic;
            const bool fixed_a = a.configuration.configs.register_type == fixed;
            const bool fixed_b = b.configuration.configs.register_type == fixed;
            if (fixed_a != fixed_b && (fixed_a ? polygon_sweep_is_separated(a, b, -normal)
                                               : polygon_sweep_is_separated(b, a, normal)))
                continue;
        }

        const auto impact_pose = [fraction](const CandidatePose& start, const CandidatePose& end) {
            return CandidatePose{start.position + (end.position - start.position) * fraction,
                                 start.rotation + (end.rotation - start.rotation) * fraction};
        };
        const auto impact_a = impact_pose(local_start_a, local_end_a);
        const auto impact_b = impact_pose(local_start_b, local_end_b);

        const WorldState::NarrowPhaseData::SweepHit sweep{toi_time,
                                                          normal,
                                                          {origin, impact_a.position,
                                                           impact_b.position, impact_a.rotation,
                                                           impact_b.rotation, true}};
        if (row)
            row->result = sweep;
        append_toi(pair, sweep);
    }

#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
    state.stats.narrow_dimension_requests += dimensions_requested;
    state.stats.narrow_dimension_loads += dimensions_loaded;
    state.stats.narrow_start_predictions += starts_predicted;
    state.stats.narrow_reused_start_predictions += predictions_reused;
#endif
    if (reuse) {
        pairs.delayed_pairs->clear();
        pairs.delayed_manifolds->clear();
        output.previous_contact_keys.clear();
    }
    previous_results.swap(next_results);
    next_results.clear();
    output.pair_result_end = state.interval.target_end;
    defer_wall_contacts(state, first_toi);
    output.prepared_time = state.schedule.time;
    output.prepared_end = state.interval.target_end;
    if (freeze_retry) {
        broad.freeze_retry = false;
        broad.reprepare_after_toi = false;
        broad.toi_changed_islands.clear();
    }
    broad.candidates.clear();

    return true;
}

bool run_narrow_phase(WorldState& state) {
    const auto& broad = state.broad_phase;
    const bool freeze_retry = broad.freeze_retry;
    const auto rigid_configs = state.body_pool.configs_column(BodyConfigs::Kind::Rigid);
    const bool cache_results =
        broad.proxy_index_ready && !freeze_retry &&
        std::any_of(rigid_configs.begin(), rigid_configs.end(), [](const BodyConfigs& config) {
            return config.register_type != BodyConfigs::RegisterType::ImmutableStatic;
        });
    return cache_results ? run_narrow_phase_candidates<true>(state)
                         : run_narrow_phase_candidates<false>(state);
}

} // namespace sqp
