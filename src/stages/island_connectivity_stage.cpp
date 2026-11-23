// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "world/world_stages.h"

#include "narrow_phase/manifold.h"

#include <algorithm>
#include <limits>
#include <numeric>
#include <variant>

namespace sqp {
namespace {

u64 body_key(BodyHandle body) {
    return (u64(body.generation) << 32) | (u64(static_cast<u8>(body.kind())) << 28) | body.index();
}

bool live_edge(const WorldState& state, const ConstraintGraph::ContactEdge& edge) {
    BodyHandle owner_a, owner_b;
    return state.body_pool.is_valid(edge.body_a) && state.body_pool.is_valid(edge.body_b) &&
           state.shape_pool.is_valid(edge.shape_a) && state.shape_pool.is_valid(edge.shape_b) &&
           state.body_pool.get_shape_owner(edge.shape_a, owner_a) &&
           state.body_pool.get_shape_owner(edge.shape_b, owner_b) && owner_a == edge.body_a &&
           owner_b == edge.body_b;
}

bool continuous(const WorldState& state, const ConstraintGraph::ContactEdge& edge) {
    if (!live_edge(state, edge))
        return false;
    BodyState a, b;
    if (!state.body_pool.get_state(edge.body_a, a) || !state.body_pool.get_state(edge.body_b, b))
        return false;
    const detail::ShapeQueryView shape_a(state.shape_pool, edge.shape_a);
    const detail::ShapeQueryView shape_b(state.shape_pool, edge.shape_b);
    ContactManifold manifold;
    if (!shape_a.is_valid() || !shape_b.is_valid() ||
        !build_contacts_from_shape_pair(shape_a, a.position, a.rotation, shape_b, b.position,
                                        b.rotation, manifold) ||
        manifold.count == 0)
        return false;
    for (i32 i = 0; i < manifold.count; ++i) {
        const Vec2 point = manifold.points[i].point;
        const Vec2 velocity_a =
            a.linear_velocity + math::cross(a.angular_velocity, point - a.position);
        const Vec2 velocity_b =
            b.linear_velocity + math::cross(b.angular_velocity, point - b.position);
        if (math::dot(velocity_a - velocity_b, manifold.normal) < config::EPSILON)
            return true;
    }
    return false;
}

u32 root(std::vector<u32>& parents, u32 index) {
    while (parents[index] != index) {
        parents[index] = parents[parents[index]];
        index = parents[index];
    }
    return index;
}

bool has_members(const Island& island) {
    if (const auto* small = std::get_if<Island::SmallIsland>(&island.bodies))
        return is_valid(small->data[0]);
    return !std::get<std::vector<BodyHandle>>(island.bodies).empty();
}

bool earlier_id(IslandId a, IslandId b) {
    return a.index != b.index ? a.index < b.index : a.generation < b.generation;
}

bool rebuild_changed_islands(WorldState& state) {
    auto& scratch = state.connectivity;
    auto sources = std::move(scratch.dirty_islands);
    for (const auto& edge : state.contact_graph.changes())
        for (BodyHandle body : {edge.body_a, edge.body_b})
            if (const auto id = state.body_pool.island_of(body))
                sources.push_back(*id);
    std::sort(sources.begin(), sources.end(), earlier_id);
    sources.erase(std::unique(sources.begin(), sources.end()), sources.end());
    std::erase_if(sources, [&](IslandId id) { return !state.island_pool.get(id); });
    std::vector<BodyHandle> bodies;
    for (IslandId id : sources) {
        const Island* island = state.island_pool.get(id);
        if (const auto* small = std::get_if<Island::SmallIsland>(&island->bodies)) {
            for (BodyHandle body : small->data) {
                if (!is_valid(body)) break;
                bodies.push_back(body);
            }
        } else {
            const auto& members = std::get<std::vector<BodyHandle>>(island->bodies);
            bodies.insert(bodies.end(), members.begin(), members.end());
        }
    }
    std::sort(bodies.begin(), bodies.end(),
              [](BodyHandle a, BodyHandle b) { return body_key(a) < body_key(b); });
    bodies.erase(std::unique(bodies.begin(), bodies.end()), bodies.end());
    const u32 count = static_cast<u32>(bodies.size());
    scratch.last_rebuilt_body_count = count;
    if (++scratch.position_epoch == 0) {
        for (auto& stamps : scratch.position_stamps)
            std::fill(stamps.begin(), stamps.end(), 0);
        scratch.position_epoch = 1;
    }
    std::vector<u32> parents(count);
    std::iota(parents.begin(), parents.end(), 0u);
    for (u32 i = 0; i < count; ++i) {
        const BodyHandle body = bodies[i];
        const auto dense = state.body_pool.dense_index(body);
        if (!dense) return stage_error(state, "rebuild_changed_islands: check failed (!dense)");
        const auto kind = static_cast<std::size_t>(body.kind());
        auto& positions = scratch.body_positions[kind];
        auto& stamps = scratch.position_stamps[kind];
        if (*dense >= positions.size()) positions.resize(*dense + 1);
        if (*dense >= stamps.size()) stamps.resize(*dense + 1, 0);
        positions[*dense] = i;
        stamps[*dense] = scratch.position_epoch;
    }
    constexpr u32 missing = std::numeric_limits<u32>::max();
    const auto position_of = [&](BodyHandle body) {
        const auto dense = state.body_pool.dense_index(body);
        if (!dense) return missing;
        const auto kind = static_cast<std::size_t>(body.kind());
        if (*dense >= scratch.position_stamps[kind].size() ||
            scratch.position_stamps[kind][*dense] != scratch.position_epoch)
            return missing;
        const u32 position = scratch.body_positions[kind][*dense];
        return position < count && bodies[position] == body ? position : missing;
    };
    for (BodyHandle body : bodies) {
        const auto edges = state.contact_graph.contacts_for(body);
        for (const auto& edge : edges) {
            const auto island_a = state.body_pool.island_of(edge.body_a);
            const auto island_b = state.body_pool.island_of(edge.body_b);
            if ((island_a && state.island_pool.frozen_until(*island_a) > state.schedule.time) ||
                (island_b && state.island_pool.frozen_until(*island_b) > state.schedule.time))
                continue;
            if (!live_edge(state, edge) || !continuous(state, edge)) {
                state.contact_graph.set_contact(
                    {edge.body_a, edge.body_b, edge.shape_a, edge.shape_b, 0, 0}, false);
                continue;
            }
            const u32 a = position_of(edge.body_a);
            const u32 b = position_of(edge.body_b);
            if (a == missing || b == missing) continue;
            const u32 ra = root(parents, a);
            const u32 rb = root(parents, b);
            if (ra != rb) parents[std::max(ra, rb)] = std::min(ra, rb);
        }
    }
    std::vector<u32> roots(count);
    std::vector<IslandId> destinations(count);
    std::vector<tick_t> times(count, 0);
    std::vector<u32> truncations(count, 0);
    std::vector<tick_t> freezes(count, 0);
    std::vector<tick_t> last_freezes(count, 0);
    for (u32 i = 0; i < count; ++i) {
        roots[i] = root(parents, i);
        const auto id = state.body_pool.island_of(bodies[i]);
        if (!id || !state.island_pool.get(*id)) return stage_error(state, "rebuild_changed_islands: check failed (!id || !state.island_pool.get(*id))");
        IslandId& destination = destinations[roots[i]];
        if (!destination.valid() || earlier_id(*id, destination)) destination = *id;
        times[roots[i]] = std::max(times[roots[i]], state.island_pool.get(*id)->time);
        truncations[roots[i]] =
            std::max(truncations[roots[i]], state.island_pool.toi_truncations(*id));
        freezes[roots[i]] =
            std::max(freezes[roots[i]], state.island_pool.frozen_until(*id));
        last_freezes[roots[i]] =
            std::max(last_freezes[roots[i]], state.island_pool.last_frozen_until(*id));
    }
    scratch.claimed_slots.clear();
    for (u32 i = 0; i < count; ++i) {
        if (roots[i] != i) continue;
        IslandId& destination = destinations[i];
        if (destination.index >= scratch.claimed_slots.size())
            scratch.claimed_slots.resize(std::size_t(destination.index) + 1, 0);
        if (scratch.claimed_slots[destination.index]) {
            destination = state.island_pool.create();
            if (!destination.valid()) return stage_error(state, "rebuild_changed_islands: check failed (!destination.valid())");
        } else {
            scratch.claimed_slots[destination.index] = 1;
        }
    }
    for (u32 i = 0; i < count; ++i) {
        const IslandId id = *state.body_pool.island_of(bodies[i]);
        const tick_t target = times[roots[i]];
        if (state.island_pool.get(id)->time < target) {
            const IslandId one[] = {id};
            if (!run_advance_islands(state, one, target)) return stage_error(state, "rebuild_changed_islands: check failed (!run_advance_islands(state, one, target))");
        }
    }
    for (u32 i = 0; i < count; ++i)
        if (roots[i] == i && !state.island_pool.set_time(destinations[i], times[i]))
            return stage_error(state, "rebuild_changed_islands: check failed (roots[i] == i && !state.island_pool.set_time(destinations[i], times[i]))");
    for (u32 i = 0; i < count; ++i) {
        const IslandId destination = destinations[roots[i]];
        if (state.body_pool.island_of(bodies[i]) != destination &&
            !state.island_pool.add_body(destination, bodies[i], state.body_pool))
            return stage_error(state, "rebuild_changed_islands: check failed (state.body_pool.island_of(bodies[i]) != destination && !state.island_pool.add_body(destination, bodies[i], state.body_pool))");
    }
    for (u32 i = 0; i < count; ++i)
        if (roots[i] == i &&
            !state.island_pool.set_toi_truncations(destinations[i], truncations[i]))
            return stage_error(state, "rebuild_changed_islands: check failed (roots[i] == i && !state.island_pool.set_toi_truncations(destinations[i], truncations[i]))");
    for (u32 i = 0; i < count; ++i)
        if (roots[i] == i &&
            !state.island_pool.set_frozen_until(destinations[i], freezes[i]))
            return stage_error(state, "rebuild_changed_islands: check failed (roots[i] == i && !state.island_pool.set_frozen_until(destinations[i], freezes[i]))");
    for (u32 i = 0; i < count; ++i)
        if (roots[i] == i &&
            !state.island_pool.set_last_frozen_until(destinations[i], last_freezes[i]))
            return stage_error(state, "rebuild_changed_islands: check failed (roots[i] == i && !state.island_pool.set_last_frozen_until(destinations[i], last_freezes[i]))");
    for (IslandId id : sources) {
        const Island* island = state.island_pool.get(id);
        if (island && !has_members(*island) &&
            (!unschedule_island(state, id) || !state.island_pool.destroy(id)))
            return stage_error(state, "rebuild_changed_islands: check failed (island && !has_members(*island) && (!unschedule_island(state, id) || !state.island_pool.destroy(id)))");
    }
    for (u32 i = 0; i < count; ++i) {
        if (roots[i] != i) continue;
        const IslandId id = destinations[i];
        if (state.island_pool.is_sleeping(id) && state.island_pool.sleep_time(id) != 0)
            continue;
        if (!wake_island(state, id)) return stage_error(state, "rebuild_changed_islands: check failed (!wake_island(state, id))");
    }
    scratch.seen_graph_version = state.contact_graph.version();
    scratch.dirty = false;
    scratch.dirty_islands.clear();
    state.contact_graph.clear_changes();
    return true;
}

} // namespace

bool run_update_contact_edges(WorldState& state) {
    for (std::size_t row = 0; row < state.pairs.pairs->size(); ++row) {
        const Pair& pair = (*state.pairs.pairs)[row];
        const bool enabled = state.pairs.manifolds->counts[row] != 0;
        state.contact_graph.set_contact(
            pair,
            enabled && continuous(state, {pair.body_a, pair.body_b, pair.shape_a, pair.shape_b}));
    }
    return true;
}

bool run_rebuild_islands(WorldState& state, bool force) {
    if (!force && !state.connectivity.dirty &&
        state.connectivity.seen_graph_version == state.contact_graph.version())
        return true;
    if (!force && !state.connectivity.full_rebuild)
        return rebuild_changed_islands(state);
    for (const auto& edge : state.contact_graph.contacts()) {
        if (!live_edge(state, edge))
            continue;
        const auto a = state.body_pool.island_of(edge.body_a);
        const auto b = state.body_pool.island_of(edge.body_b);
        if ((a && state.island_pool.frozen_until(*a) > state.schedule.time) ||
            (b && state.island_pool.frozen_until(*b) > state.schedule.time))
            continue;
        if (!a || !b || *a == *b)
            continue;
        const Island* island_a = state.island_pool.get(*a);
        const Island* island_b = state.island_pool.get(*b);
        if (!island_a || !island_b)
            return stage_error(state, "run_rebuild_islands: check failed (!island_a || !island_b)");
        if (island_a->time < island_b->time) {
            const IslandId one[] = {*a};
            if (!run_advance_islands(state, one, island_b->time))
                return stage_error(state, "run_rebuild_islands: check failed (!run_advance_islands(state, one, island_b->time))");
        } else if (island_b->time < island_a->time) {
            const IslandId one[] = {*b};
            if (!run_advance_islands(state, one, island_a->time))
                return stage_error(state, "run_rebuild_islands: check failed (!run_advance_islands(state, one, island_a->time))");
        }
    }
    state.contact_graph.erase_if([&](const auto& edge) { return !continuous(state, edge); });
    auto& scratch = state.connectivity;
    state.body_pool.get_handles(scratch.bodies);
    std::erase_if(scratch.bodies, [&](BodyHandle body) {
        BodyConfigs::RegisterType type;
        return !state.body_pool.get_register_type(body, type) ||
               type == BodyConfigs::RegisterType::ImmutableStatic;
    });
    std::sort(scratch.bodies.begin(), scratch.bodies.end(),
              [](BodyHandle a, BodyHandle b) { return body_key(a) < body_key(b); });
    const u32 count = static_cast<u32>(scratch.bodies.size());
    scratch.last_rebuilt_body_count = count;
    scratch.parents.resize(count);
    std::iota(scratch.parents.begin(), scratch.parents.end(), 0u);
    constexpr u32 missing = std::numeric_limits<u32>::max();
    for (BodyConfigs::Kind kind :
         {BodyConfigs::Kind::Rigid, BodyConfigs::Kind::LiteCircle, BodyConfigs::Kind::AabbBox})
        scratch.body_positions[static_cast<std::size_t>(kind)].assign(
            state.body_pool.configs_column(kind).size(), missing);
    for (u32 i = 0; i < count; ++i) {
        const BodyHandle body = scratch.bodies[i];
        const auto dense = state.body_pool.dense_index(body);
        auto& positions = scratch.body_positions[static_cast<std::size_t>(body.kind())];
        if (!dense || *dense >= positions.size())
            return stage_error(state, "run_rebuild_islands: check failed (!dense || *dense >= positions.size())");
        positions[*dense] = i;
    }
    const auto position_of = [&](BodyHandle body) {
        const auto dense = state.body_pool.dense_index(body);
        const auto kind = static_cast<std::size_t>(body.kind());
        if (!dense || kind >= scratch.body_positions.size() ||
            *dense >= scratch.body_positions[kind].size())
            return missing;
        return scratch.body_positions[kind][*dense];
    };
    for (const auto& edge : state.contact_graph.contacts()) {
        const u32 a = position_of(edge.body_a);
        const u32 b = position_of(edge.body_b);
        const auto island_a = state.body_pool.island_of(edge.body_a);
        const auto island_b = state.body_pool.island_of(edge.body_b);
        if ((island_a && state.island_pool.frozen_until(*island_a) > state.schedule.time) ||
            (island_b && state.island_pool.frozen_until(*island_b) > state.schedule.time))
            continue;
        if (a == missing || b == missing)
            continue;
        u32 ra = root(scratch.parents, a);
        u32 rb = root(scratch.parents, b);
        if (ra != rb)
            scratch.parents[std::max(ra, rb)] = std::min(ra, rb);
    }
    scratch.roots.resize(count);
    for (u32 i = 0; i < count; ++i)
        scratch.roots[i] = root(scratch.parents, i);
    state.island_pool.get_ids(scratch.old_islands);
    scratch.destinations.assign(count, {});
    scratch.component_toi_truncations.assign(count, 0);
    std::vector<tick_t> component_freezes(count, 0);
    std::vector<tick_t> component_last_freezes(count, 0);
    scratch.claimed_slots.assign(state.schedule.entries.size(), 0);
    for (u32 i = 0; i < count; ++i) {
        const auto id = state.body_pool.island_of(scratch.bodies[i]);
        if (!id || !state.island_pool.get(*id))
            return stage_error(state, "run_rebuild_islands: check failed (!id || !state.island_pool.get(*id))");
        const u32 component = scratch.roots[i];
        const u32 source_truncations = state.island_pool.toi_truncations(*id);
        scratch.component_toi_truncations[component] =
            std::max(scratch.component_toi_truncations[component], source_truncations);
        component_freezes[component] =
            std::max(component_freezes[component], state.island_pool.frozen_until(*id));
        component_last_freezes[component] =
            std::max(component_last_freezes[component],
                     state.island_pool.last_frozen_until(*id));
        if (id->index >= scratch.claimed_slots.size())
            scratch.claimed_slots.resize(std::size_t(id->index) + 1, 0);
        if (!scratch.destinations[component].valid() ||
            earlier_id(*id, scratch.destinations[component]))
            scratch.destinations[component] = *id;
    }
    // One old ID can survive only one component. Components are visited in body-handle order.
    for (u32 i = 0; i < count; ++i) {
        if (scratch.roots[i] != i)
            continue;
        IslandId& destination = scratch.destinations[i];
        if (scratch.claimed_slots[destination.index]) {
            destination = state.island_pool.create();
            if (!destination.valid())
                return stage_error(state, "run_rebuild_islands: check failed (!destination.valid())");
        } else {
            scratch.claimed_slots[destination.index] = 1;
        }
    }
    // Synchronize all source islands in each component before moving their members.
    for (u32 i = 0; i < count; ++i) {
        if (scratch.roots[i] != i)
            continue;
        tick_t target = 0;
        for (u32 j = 0; j < count; ++j)
            if (scratch.roots[j] == i) {
                const auto id = state.body_pool.island_of(scratch.bodies[j]);
                target = std::max(target, state.island_pool.get(*id)->time);
            }
        for (u32 j = 0; j < count; ++j)
            if (scratch.roots[j] == i) {
                const auto id = state.body_pool.island_of(scratch.bodies[j]);
                if (state.island_pool.get(*id)->time < target) {
                    const IslandId one[] = {*id};
                    if (!run_advance_islands(state, one, target))
                        return stage_error(state, "run_rebuild_islands: check failed (!run_advance_islands(state, one, target))");
                }
            }
        if (!state.island_pool.set_time(scratch.destinations[i], target))
            return stage_error(state, "run_rebuild_islands: check failed (!state.island_pool.set_time(scratch.destinations[i], target))");
    }
    for (u32 i = 0; i < count; ++i) {
        const IslandId destination = scratch.destinations[scratch.roots[i]];
        if (state.body_pool.island_of(scratch.bodies[i]) != destination &&
            !state.island_pool.add_body(destination, scratch.bodies[i], state.body_pool))
            return stage_error(state, "run_rebuild_islands: check failed (state.body_pool.island_of(scratch.bodies[i]) != destination && !state.island_pool.add_body(destination, scratch.bodies[i], state.body_pool))");
    }
    for (u32 i = 0; i < count; ++i)
        if (scratch.roots[i] == i &&
            !state.island_pool.set_toi_truncations(scratch.destinations[i],
                                                   scratch.component_toi_truncations[i]))
            return stage_error(state, "run_rebuild_islands: check failed (scratch.roots[i] == i && !state.island_pool.set_toi_truncations(scratch.destinations[i], scratch.component_toi_truncations[i]))");
    for (u32 i = 0; i < count; ++i)
        if (scratch.roots[i] == i &&
            !state.island_pool.set_frozen_until(scratch.destinations[i], component_freezes[i]))
            return stage_error(state, "run_rebuild_islands: check failed (scratch.roots[i] == i && !state.island_pool.set_frozen_until(scratch.destinations[i], component_freezes[i]))");
    for (u32 i = 0; i < count; ++i)
        if (scratch.roots[i] == i &&
            !state.island_pool.set_last_frozen_until(scratch.destinations[i],
                                                    component_last_freezes[i]))
            return stage_error(state, "run_rebuild_islands: check failed (scratch.roots[i] == i && !state.island_pool.set_last_frozen_until(scratch.destinations[i], component_last_freezes[i]))");
    for (IslandId id : scratch.old_islands) {
        const Island* island = state.island_pool.get(id);
        if (island && !has_members(*island))
            if (!unschedule_island(state, id) || !state.island_pool.destroy(id))
                return stage_error(state, "run_rebuild_islands: check failed (!unschedule_island(state, id) || !state.island_pool.destroy(id))");
    }
    for (u32 i = 0; i < count; ++i) {
        if (scratch.roots[i] != i)
            continue;
        const IslandId id = scratch.destinations[i];
        u8 bucket;
        std::size_t position;
        if (!state.island_pool.get_schedule_location(id, bucket, position))
            return stage_error(state, "run_rebuild_islands: check failed (!state.island_pool.get_schedule_location(id, bucket, position))");
        if (bucket == 0 && state.island_pool.sleep_time(id) != 0)
            continue;
        if (!wake_island(state, id))
            return stage_error(state, "run_rebuild_islands: check failed (!wake_island(state, id))");
    }
    scratch.seen_graph_version = state.contact_graph.version();
    scratch.dirty = false;
    scratch.full_rebuild = false;
    scratch.dirty_islands.clear();
    state.contact_graph.clear_changes();
    return true;
}

} // namespace sqp
