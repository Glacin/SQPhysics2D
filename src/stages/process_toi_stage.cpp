// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "world/world_stages.h"

#include "narrow_phase/manifold.h"
#include "solver/integration.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <variant>

namespace sqp {
namespace {

using ToiEvent = WorldState::ProcessToiData::ToiEvent;
using ToiEventCompare = WorldState::ProcessToiData::ToiEventCompare;

bool current_version(const WorldState& state, BodyHandle body, u64 expected) {
    const auto index = state.body_pool.dense_index(body);
    if (!index)
        return false;
    const auto versions = state.body_pool.motion_version_column(body.kind());
    return *index < versions.size() && versions[*index] == expected;
}

bool event_is_current(const WorldState& state, const ToiEvent& event) {
    BodyHandle owner_a, owner_b;
    return current_version(state, event.pair.body_a, event.pair.motion_version_a) &&
           current_version(state, event.pair.body_b, event.pair.motion_version_b) &&
           state.shape_pool.is_valid(event.pair.shape_a) &&
           state.shape_pool.is_valid(event.pair.shape_b) &&
           state.body_pool.get_shape_owner(event.pair.shape_a, owner_a) &&
           state.body_pool.get_shape_owner(event.pair.shape_b, owner_b) &&
           owner_a == event.pair.body_a && owner_b == event.pair.body_b;
}

bool event_is_after_islands(const WorldState& state, const ToiEvent& event) {
    const auto id_a = state.body_pool.island_of(event.pair.body_a);
    const auto id_b = state.body_pool.island_of(event.pair.body_b);
    const Island* island_a = id_a ? state.island_pool.get(*id_a) : nullptr;
    const Island* island_b = id_b ? state.island_pool.get(*id_b) : nullptr;
    BodyConfigs::RegisterType type_a, type_b;
    if (!state.body_pool.get_register_type(event.pair.body_a, type_a) ||
        !state.body_pool.get_register_type(event.pair.body_b, type_b))
        return false;
    const bool frozen_a = id_a && state.island_pool.frozen_until(*id_a) > state.schedule.time;
    const bool frozen_b = id_b && state.island_pool.frozen_until(*id_b) > state.schedule.time;
    return (island_a ? (frozen_a ? event.time >= state.schedule.time &&
                                      event.time < state.island_pool.frozen_until(*id_a)
                                : event.time > island_a->time)
                     : type_a == BodyConfigs::RegisterType::ImmutableStatic) &&
           (island_b ? (frozen_b ? event.time >= state.schedule.time &&
                                      event.time < state.island_pool.frozen_until(*id_b)
                                : event.time > island_b->time)
                     : type_b == BodyConfigs::RegisterType::ImmutableStatic);
}

bool propagate_frozen_contact(WorldState& state, const Pair& pair, tick_t time,
                              bool& propagated) {
    propagated = false;
    const auto a = state.body_pool.island_of(pair.body_a);
    const auto b = state.body_pool.island_of(pair.body_b);
    const bool frozen_a = a && state.island_pool.frozen_until(*a) > time;
    const bool frozen_b = b && state.island_pool.frozen_until(*b) > time;
    if (!frozen_a && !frozen_b)
        return true;
    if (frozen_a && b && !frozen_b && !freeze_island_until_next_tick(state, *b, time))
        return stage_error(state, "propagate_frozen_contact: check failed (frozen_a && b && !frozen_b && !freeze_island_until_next_tick(state, *b, time))");
    if (frozen_b && a && !frozen_a && !freeze_island_until_next_tick(state, *a, time))
        return stage_error(state, "propagate_frozen_contact: check failed (frozen_b && a && !frozen_a && !freeze_island_until_next_tick(state, *a, time))");
    propagated = true;
    return true;
}

bool event_has_budget(const WorldState& state, const ToiEvent& event) {
    const u32 limit = state.settings.interval.max_truncations_per_tick;
    for (BodyHandle body : {event.pair.body_a, event.pair.body_b}) {
        BodyConfigs::RegisterType type;
        if (!state.body_pool.get_register_type(body, type))
            return false;
        if (type == BodyConfigs::RegisterType::ImmutableStatic)
            continue;
        const auto island = state.body_pool.island_of(body);
        if (!island || !state.island_pool.get(*island) ||
            state.island_pool.toi_truncations(*island) >= limit)
            return false;
    }
    return true;
}

bool pose_at(const WorldState& state, BodyHandle body, tick_t time, Vec2& position,
             real& rotation) {
    BodyState body_state;
    BodyConfigs::RegisterType registration;
    const auto id = state.body_pool.island_of(body);
    const Island* island = id ? state.island_pool.get(*id) : nullptr;
    if (!state.body_pool.get_state(body, body_state) ||
        !state.body_pool.get_register_type(body, registration) ||
        (registration != BodyConfigs::RegisterType::ImmutableStatic && !island))
        return false;
    position = body_state.position;
    rotation = body_state.rotation;
    if (registration != BodyConfigs::RegisterType::ImmutableStatic &&
        state.island_pool.frozen_until(*id) <= time) {
        if (island->time >= time)
            return false;
        predict_body_state(body_state.position, body_state.rotation, body_state.linear_velocity,
                           body_state.angular_velocity, {}, real(0), time - island->time, position,
                           rotation);
    }
    return true;
}

struct SupportFeature {
    Vec2 low, high;
    ContactFeatureRef feature;
};

bool support_feature(const detail::ShapeQueryView& shape, Vec2 position, real rotation,
                     Vec2 direction, Vec2 tangent, real tolerance, SupportFeature& out) {
    if (!shape.is_valid())
        return false;
    if (const auto* circle = std::get_if<Circle>(&shape.geometry())) {
        if (circle->radius <= real(0))
            return false;
        out.low = out.high = position + math::rotate_vector(shape.local_center, rotation) +
                             direction * circle->radius;
        out.feature.type = ContactFeatureType::Circle;
        return true;
    }
    const auto* polygon = std::get_if<ConvexPoly>(&shape.geometry());
    if (!polygon || polygon->local.size() < 3)
        return false;
    const Vec2 center = position + math::rotate_vector(polygon->center_offset(), rotation);
    const auto world_vertex = [&](std::size_t index) {
        return center + math::rotate_vector(polygon->local[index], rotation);
    };
    real maximum = -std::numeric_limits<real>::infinity();
    u16 best_index = 0;
    for (std::size_t i = 0; i < polygon->local.size(); ++i) {
        const real projection = math::dot(world_vertex(i), direction);
        if (projection > maximum) {
            maximum = projection;
            best_index = static_cast<u16>(i);
        }
    }
    if (!std::isfinite(maximum))
        return false;
    for (std::size_t i = 0; i < polygon->local.size(); ++i) {
        const Vec2 first = world_vertex(i);
        const Vec2 second = world_vertex((i + 1) % polygon->local.size());
        const real edge_tolerance = std::max(tolerance, (second - first).length() * real(1e-4));
        if (maximum - math::dot(first, direction) > edge_tolerance ||
            maximum - math::dot(second, direction) > edge_tolerance)
            continue;
        out.low = first;
        out.high = second;
        if (math::dot(out.low, tangent) > math::dot(out.high, tangent))
            std::swap(out.low, out.high);
        out.feature.type = ContactFeatureType::Edge;
        out.feature.index = static_cast<u16>(i);
        return true;
    }
    out.low = out.high = world_vertex(best_index);
    out.feature.type = ContactFeatureType::Vertex;
    out.feature.index = best_index;
    return true;
}

Vec2 at_tangent(const SupportFeature& feature, Vec2 tangent, real coordinate) {
    const real low = math::dot(feature.low, tangent);
    const real high = math::dot(feature.high, tangent);
    if (high <= low)
        return feature.low;
    const real fraction = std::clamp((coordinate - low) / (high - low), real(0), real(1));
    return feature.low + (feature.high - feature.low) * fraction;
}

bool zero_depth_contact(const detail::ShapeQueryView& shape_a, Vec2 position_a, real rotation_a,
                        const detail::ShapeQueryView& shape_b, Vec2 position_b, real rotation_b,
                        const ToiEvent& event, real tolerance, ContactManifold& manifold) {
    const real normal_length_squared = event.normal.length_squared();
    if (!std::isfinite(normal_length_squared) ||
        normal_length_squared <= config::EPSILON * config::EPSILON)
        return false;
    const Vec2 normal = -event.normal / math::sqrt(normal_length_squared);
    const Vec2 tangent = math::perp_ccw(normal);
    SupportFeature a, b;
    const real feature_tolerance = config::EPSILON * real(8);
    if (!support_feature(shape_a, position_a, rotation_a, normal, tangent, feature_tolerance, a) ||
        !support_feature(shape_b, position_b, rotation_b, -normal, tangent, feature_tolerance, b))
        return false;
    const real low = std::max(math::dot(a.low, tangent), math::dot(b.low, tangent));
    const real high = std::min(math::dot(a.high, tangent), math::dot(b.high, tangent));
    if (low > high + tolerance)
        return false;
    const real coordinate = (low + high) * real(0.5);
    const Vec2 point_a = at_tangent(a, tangent, coordinate);
    const Vec2 point_b = at_tangent(b, tangent, coordinate);
    const real separation = math::dot(point_b - point_a, normal);
    if (!std::isfinite(separation) || math::abs(separation) > tolerance)
        return false;
    manifold = {};
    manifold.count = 1;
    manifold.normal = normal;
    manifold.points[0].point = (point_a + point_b) * real(0.5);
    manifold.points[0].key.feature_a = a.feature;
    manifold.points[0].key.feature_b = b.feature;
    return true;
}

bool event_manifold(const WorldState& state, const ToiEvent& event, ContactManifold& manifold,
                    WorldState::ToiLocalContact& local) {
    BodyHandle owner_a, owner_b;
    if (!state.shape_pool.is_valid(event.pair.shape_a) ||
        !state.shape_pool.is_valid(event.pair.shape_b) ||
        !state.body_pool.get_shape_owner(event.pair.shape_a, owner_a) ||
        !state.body_pool.get_shape_owner(event.pair.shape_b, owner_b) ||
        owner_a != event.pair.body_a || owner_b != event.pair.body_b)
        return false;
    Vec2 position_a, position_b;
    real rotation_a = 0, rotation_b = 0;
    if (!pose_at(state, event.pair.body_a, event.time, position_a, rotation_a) ||
        !pose_at(state, event.pair.body_b, event.time, position_b, rotation_b))
        return false;
    const detail::ShapeQueryView shape_a(state.shape_pool, event.pair.shape_a);
    const detail::ShapeQueryView shape_b(state.shape_pool, event.pair.shape_b);
    const auto& impact = event.impact;
    local.origin = impact.valid ? impact.origin : position_a;
    const Vec2 local_a = impact.valid ? impact.position_a : Vec2{};
    const Vec2 local_b = impact.valid ? impact.position_b : position_b - position_a;
    const real angle_a = impact.valid ? impact.rotation_a : rotation_a;
    const real angle_b = impact.valid ? impact.rotation_b : rotation_b;
    if (!build_contacts_from_shape_pair(shape_a, local_a, angle_a, shape_b, local_b, angle_b,
                                        local.manifold)) {
        // This works, but the real reason is not finded yet.
        // Manifold is discarded in some situations if I set it to config::EPSILON
        constexpr real tolerance = real(1e-4);
        if (!zero_depth_contact(shape_a, local_a, angle_a, shape_b, local_b, angle_b, event,
                                tolerance, local.manifold))
            return false;
    }
    local.valid = true;
    manifold = local.manifold;
    for (i32 i = 0; i < manifold.count; ++i)
        manifold.points[i].point += local.origin;
    return true;
}

} // namespace

// Might need to refactor the whole stage, TOI=0 is really complex when I'am trying to step to exact
// time. it might freeze the body.
bool run_process_toi(WorldState& state, tick_t bound, bool& selected_time) {
    selected_time = false;
    auto& input = state.narrow_phase;
    auto& output = state.process_toi;
    auto& pairs = state.pairs;
    auto& interval = state.interval;
    if (pairs.pairs->size() != pairs.manifolds->size() ||
        (!pairs.toi_local_contacts.empty() &&
         pairs.toi_local_contacts.size() != pairs.pairs->size()) ||
        pairs.delayed_pairs->size() != pairs.delayed_manifolds->size() ||
        (!pairs.delayed_toi_local_contacts.empty() &&
         pairs.delayed_toi_local_contacts.size() != pairs.delayed_pairs->size()) ||
        (!interval.microstep &&
         (!pairs.delayed_pairs->empty() || !pairs.delayed_manifolds->empty() ||
          !pairs.delayed_toi_local_contacts.empty() || !interval.delayed_islands->empty())) ||
        input.toi_pairs.size() != input.toi_times.size() ||
        input.toi_pairs.size() != input.toi_normals.size() ||
        (!input.toi_impacts.empty() && input.toi_impacts.size() != input.toi_pairs.size()) ||
        input.toi_pairs.size() > std::numeric_limits<u64>::max() - output.next_sequence)
        return stage_error(state, "run_process_toi: TOI or contact pair columns have inconsistent sizes");

    const ToiEventCompare compare;
    const auto swap_active_delayed = [&] {
        pairs.pairs.swap(pairs.delayed_pairs);
        pairs.manifolds.swap(pairs.delayed_manifolds);
        pairs.toi_local_contacts.swap(pairs.delayed_toi_local_contacts);
        std::swap(pairs.original_manifolds, pairs.delayed_original_manifolds);
        std::swap(pairs.original_manifolds_valid, pairs.delayed_original_manifolds_valid);
        interval.islands.swap(interval.delayed_islands);
    };
    for (std::size_t i = 0; i < input.toi_pairs.size(); ++i) {
        const bool duplicate =
            std::any_of(output.toi_heap.begin(), output.toi_heap.end(), [&](const ToiEvent& event) {
                return event.time == input.toi_times[i] && event.pair == input.toi_pairs[i];
            });
        if (!duplicate) {
            output.toi_heap.push_back(
                {input.toi_times[i], input.toi_pairs[i], input.toi_normals[i],
                 output.next_sequence++,
                 input.toi_impacts.empty() ? WorldState::ToiLocalImpact{} : input.toi_impacts[i]});
            std::push_heap(output.toi_heap.begin(), output.toi_heap.end(), compare);
        }
    }
    input.toi_pairs.clear();
    input.toi_times.clear();
    input.toi_normals.clear();
    input.toi_impacts.clear();

    if (!interval.microstep) {
        for (std::size_t i = 0; i < pairs.pairs->size(); ++i) {
            if (pairs.manifolds->counts[i] == 0)
                continue;
            bool propagated = false;
            if (!propagate_frozen_contact(state, (*pairs.pairs)[i], interval.target_end,
                                          propagated))
                return stage_error(state, "run_process_toi: check failed (!propagate_frozen_contact(state, (*pairs.pairs)[i], interval.target_end, propagated))");
            if (propagated) {
                std::erase_if(*interval.islands, [&](IslandId id) {
                    return state.island_pool.frozen_until(id) > state.schedule.time &&
                           state.island_pool.get(id)->time > interval.target_end;
                });
                state.schedule.prepared = false;
                return true;
            }
        }
    }

    // Preserve a microstep that was already selected by an earlier call.
    if (state.interval.microstep)
        return true;

    while (!output.toi_heap.empty() && output.toi_heap.front().time <= bound) {
        const tick_t earliest = output.toi_heap.front().time;
        bool valid_time = false;
        const auto restore_ordinary = [&] {
            if (!valid_time)
                return;
            pairs.pairs->clear();
            pairs.manifolds->clear();
            pairs.toi_local_contacts.clear();
            pairs.original_manifolds.clear();
            pairs.original_manifolds_valid = false;
            swap_active_delayed();
        };
        do {
            std::pop_heap(output.toi_heap.begin(), output.toi_heap.end(), compare);
            ToiEvent event = output.toi_heap.back();
            output.toi_heap.pop_back();
            if (event.time < state.schedule.time || !event_is_current(state, event))
                continue;
            if (!event_is_after_islands(state, event))
                continue;
            ContactManifold manifold;
            WorldState::ToiLocalContact local;
            if (!event_manifold(state, event, manifold, local))
                continue;
            bool propagated = false;
            if (!propagate_frozen_contact(state, event.pair, event.time, propagated))
                return stage_error(state, "run_process_toi: check failed (!propagate_frozen_contact(state, event.pair, event.time, propagated))");
            if (propagated) {
                restore_ordinary();
                std::erase_if(*interval.islands, [&](IslandId id) {
                    return state.island_pool.frozen_until(id) > state.schedule.time &&
                           state.island_pool.get(id)->time > interval.target_end;
                });
                state.schedule.prepared = false;
                return true;
            }
            if (!event_has_budget(state, event))
                continue;
            if (!valid_time) {
                swap_active_delayed();
            }
            valid_time = true;
            if (!pairs.manifolds->append(manifold)) {
                restore_ordinary();
                return stage_error(state, "run_process_toi: check failed (!pairs.manifolds->append(manifold))");
            }
            pairs.pairs->push_back(event.pair);
            pairs.toi_local_contacts.push_back(local);
        } while (!output.toi_heap.empty() && output.toi_heap.front().time == earliest);

        if (valid_time) {
            state.interval.microstep = true;
            state.interval.tick = false;
            state.interval.microstep_time = earliest;
            selected_time = true;
            return true;
        }
    }
    state.interval.microstep = false;
    return true;
}

} // namespace sqp
