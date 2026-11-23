// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "world/world_stages.h"

#include "solver/integration.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <variant>

namespace sqp {
namespace {

bool valid_layout(const WorldState::PairsData& pairs) {
    if (!pairs.pairs || !pairs.manifolds)
        return false;
    const std::size_t count = pairs.pairs->size();
    if (count > std::numeric_limits<std::size_t>::max() / 2)
        return false;
    const auto& m = *pairs.manifolds;
    const std::size_t points = count * 2;
    return m.counts.size() == count && m.normals_x.size() == count && m.normals_y.size() == count &&
           m.keys.size() == points && m.points_x.size() == points && m.points_y.size() == points &&
           m.depths.size() == points && pairs.combined_restitutions.size() == count &&
           pairs.combined_dynamic_frictions.size() == count &&
           pairs.combined_static_frictions.size() == count &&
           pairs.combined_damping_ratios.size() == count &&
           pairs.combined_contact_frequencies.size() == count &&
           pairs.combined_max_penetration_correction_speeds.size() == count;
}

#ifdef SQP_DEBUG_VALIDATION
bool valid_parameters(const WorldState::PairsData& pairs, std::size_t i) {
    const real restitution = pairs.combined_restitutions[i];
    const real dynamic_friction = pairs.combined_dynamic_frictions[i];
    const real static_friction = pairs.combined_static_frictions[i];
    return std::isfinite(restitution) && restitution >= real(0) && restitution <= real(1) &&
           std::isfinite(dynamic_friction) && dynamic_friction >= real(0) &&
           std::isfinite(static_friction) && static_friction >= dynamic_friction &&
           std::isfinite(pairs.combined_damping_ratios[i]) &&
           pairs.combined_damping_ratios[i] >= real(0) &&
           std::isfinite(pairs.combined_contact_frequencies[i]) &&
           pairs.combined_contact_frequencies[i] >= real(0) &&
           std::isfinite(pairs.combined_max_penetration_correction_speeds[i]) &&
           pairs.combined_max_penetration_correction_speeds[i] >= real(0);
}
#endif

bool valid_manifold(const ContactManifold& manifold) {
    if (manifold.count < 0 || manifold.count > 2)
        return false;
#ifdef SQP_DEBUG_VALIDATION
    if (manifold.count == 0)
        return true;
    if (!manifold.normal.isfinite() ||
        !(manifold.normal.length_squared() > config::EPSILON * config::EPSILON))
        return false;
    for (i32 point = 0; point < manifold.count; ++point)
        if (!manifold.points[point].point.isfinite() ||
            !std::isfinite(manifold.points[point].depth))
            return false;
#endif
    return true;
}

bool same_contact_pair(const WorldState::SolveData::ContactCacheEntry& entry,
                       const Pair& pair) {
    return entry.body_a == pair.body_a && entry.body_b == pair.body_b &&
           entry.shape_a == pair.shape_a && entry.shape_b == pair.shape_b;
}

void restore_contact_impulses(const WorldState::SolveData& solve, const Pair& pair,
                              const ContactManifold& manifold, tick_t start, real step_seconds,
                              ContactConstraint& constraint) {
    const auto found = std::find_if(solve.contact_cache.begin(), solve.contact_cache.end(),
                                    [&](const auto& entry) {
                                        return same_contact_pair(entry, pair);
                                    });
    if (found == solve.contact_cache.end() || found->time != start ||
        math::dot(found->normal, constraint.normal) < real(0.95))
        return;
    const real ratio = found->step_seconds > real(0)
                           ? step_seconds / found->step_seconds
                           : real(0);
    const real warm_scale = std::isfinite(ratio) && ratio > real(0) ? ratio : real(0);
    bool used[2] = {};
    for (i32 i = 0; i < constraint.count; ++i) {
        for (u8 j = 0; j < found->count; ++j) {
            if (used[j] || found->keys[j] != manifold.points[i].key)
                continue;
            auto& point = constraint.points[i];
            point.normal_impulse = std::max(real(0), found->normal_impulses[j] * warm_scale);
            point.tangent_impulse = std::clamp(
                found->tangent_impulses[j] * warm_scale,
                -constraint.static_friction * point.normal_impulse,
                constraint.static_friction * point.normal_impulse);
            point.restitution_velocity = 0;
            used[j] = true;
            break;
        }
    }
}

void save_contact_impulses(WorldState::SolveData& solve, const Pair& pair,
                          const ContactManifold& manifold, const ContactConstraint& constraint,
                          tick_t target, real step_seconds) {
    auto found = std::find_if(solve.contact_cache.begin(), solve.contact_cache.end(),
                              [&](const auto& entry) {
                                  return same_contact_pair(entry, pair);
                              });
    WorldState::SolveData::ContactCacheEntry entry;
    entry.body_a = pair.body_a;
    entry.body_b = pair.body_b;
    entry.shape_a = pair.shape_a;
    entry.shape_b = pair.shape_b;
    entry.normal = constraint.normal;
    entry.time = target;
    entry.step_seconds = step_seconds;
    bool supported = false;
    for (i32 i = 0; i < constraint.count; ++i)
        supported |= constraint.points[i].peak_normal_impulse > real(0);
    for (i32 i = 0; i < constraint.count; ++i) {
        const auto& point = constraint.points[i];
        if (!supported || point.current_separation > real(0))
            continue;
        const u8 slot = entry.count++;
        entry.keys[slot] = manifold.points[i].key;
        entry.normal_impulses[slot] = point.normal_impulse;
        entry.tangent_impulses[slot] = point.tangent_impulse;
    }
    if (found != solve.contact_cache.end()) {
        if (entry.count == 0)
            solve.contact_cache.erase(found);
        else
            *found = entry;
    } else if (entry.count > 0) {
        solve.contact_cache.push_back(entry);
    }
}

void prune_contact_cache(WorldState& state) {
    auto& cache = state.solve.contact_cache;
    std::erase_if(cache, [&](const auto& entry) {
        if (!state.body_pool.is_valid(entry.body_a) ||
            !state.body_pool.is_valid(entry.body_b) ||
            !state.shape_pool.is_valid(entry.shape_a) ||
            !state.shape_pool.is_valid(entry.shape_b))
            return true;
        const auto advanced = [&](BodyHandle body) {
            const auto island = state.body_pool.island_of(body);
            const auto* value = island ? state.island_pool.get(*island) : nullptr;
            return value && value->time > entry.time;
        };
        return advanced(entry.body_a) || advanced(entry.body_b);
    });
}

bool add_island(WorldState& state, IslandId id, tick_t target) {
    const Island* island = state.island_pool.get(id);
    if (!island || island->time > target)
        return stage_error(state, "add_island: check failed (!island || island->time > target)");
    auto& solve = state.solve;
    if (id.index >= solve.island_stamps.size())
        solve.island_stamps.resize(std::size_t(id.index) + 1, 0);
    if (solve.island_stamps[id.index] != solve.selection_epoch) {
        solve.island_stamps[id.index] = solve.selection_epoch;
        solve.islands.push_back(id);
    }
    return true;
}

bool collect_bodies(WorldState& state, tick_t target) {
    auto& solve = state.solve;
    for (IslandId id : solve.islands) {
        const Island* island = state.island_pool.get(id);
        if (!island || island->time > target)
            return stage_error(state, "collect_bodies: check failed (!island || island->time > target)");
        const auto append = [&](BodyHandle body) {
            if (!state.body_pool.is_valid(body) || state.body_pool.island_of(body) != id)
                return stage_error(state, "collect_bodies: check failed (!state.body_pool.is_valid(body) || state.body_pool.island_of(body) != id)");
            const auto dense = state.body_pool.dense_index(body);
            const auto kind = static_cast<std::size_t>(body.kind());
            if (!dense || kind >= solve.body_stamps.size())
                return stage_error(state, "collect_bodies: check failed (!dense || kind >= solve.body_stamps.size())");
            auto& stamps = solve.body_stamps[kind];
            if (*dense >= stamps.size())
                stamps.resize(std::size_t(*dense) + 1, 0);
            if (stamps[*dense] == solve.selection_epoch) {
#ifdef SQP_DEBUG_VALIDATION
                return stage_error(state, "collect_bodies: check failed (stamps[*dense] == solve.selection_epoch)");
#else
                return true;
#endif
            }
            stamps[*dense] = solve.selection_epoch;
            solve.bodies.push_back(body);
            solve.body_start_times.push_back(island->time);
            return true;
        };
        if (const auto* small = std::get_if<Island::SmallIsland>(&island->bodies)) {
            for (BodyHandle body : small->data) {
                if (!is_valid(body))
                    break;
                if (!append(body))
                    return stage_error(state, "collect_bodies: check failed (!append(body))");
            }
        } else {
            for (BodyHandle body : std::get<std::vector<BodyHandle>>(island->bodies))
                if (!append(body))
                    return stage_error(state, "collect_bodies: check failed (!append(body))");
        }
    }
    return true;
}

bool valid_body(const WorldState& state, BodyHandle body, tick_t start, tick_t target) {
    if (start > target)
        return false;
#ifdef SQP_DEBUG_VALIDATION
    BodyState value;
    BodyConfiguration config;
    BodyProperties properties;
    Vec2 force;
    real torque;
    if (!state.body_pool.get_state(body, value) || !value.is_alive() ||
        !state.body_pool.get_configuration(body, config) ||
        !state.body_pool.get_properties(body, properties) ||
        !state.body_pool.get_constant_force(body, force, torque) || !value.position.isfinite() ||
        !value.linear_velocity.isfinite() || !std::isfinite(value.rotation) ||
        !std::isfinite(value.angular_velocity) || !force.isfinite() || !std::isfinite(torque) ||
        !std::isfinite(properties.inv_mass) || !std::isfinite(properties.inv_inertia))
        return false;
#else
    (void)state;
    (void)body;
#endif
    return std::isfinite(math::tu_to_seconds<real>(target - start));
}

void bump_motion_version(WorldState& state, BodyHandle body) {
    const auto dense = state.body_pool.dense_index(body);
    if (!dense)
        return;
    auto versions = state.body_pool.motion_version_column(body.kind());
    if (*dense < versions.size())
        ++versions[*dense];
}

bool advance_velocity(WorldState& state, BodyHandle body, real seconds) {
    BodyState value;
    BodyConfiguration config;
    BodyProperties properties;
    Vec2 force;
    real torque;
    if (!state.body_pool.get_state(body, value) ||
        !state.body_pool.get_configuration(body, config) ||
        !state.body_pool.get_properties(body, properties) ||
        !state.body_pool.get_constant_force(body, force, torque))
        return stage_error(state, "advance_velocity: check failed (!state.body_pool.get_state(body, value) || !state.body_pool.get_configuration(body, config) || !state.body_pool.get_properties(body, properties) || !state.body_pool.get_constant_force(body, force, torque))");
    if (config.configs.register_type == BodyConfigs::RegisterType::Dynamic) {
        const Vec2 gravity =
            gravity_with_body_overrides(state.settings.gravity, config.configs.override_flags);
        integrate_body_velocity_seconds(value.linear_velocity, value.angular_velocity,
                                        gravity + force * properties.inv_mass,
                                        torque * properties.inv_inertia, seconds);
        const auto flags = config.configs.override_flags;
        if (math::has_any_flags(flags, BodyConfigs::OverrideFlags::LockLinearX))
            value.linear_velocity.x = real(0);
        if (math::has_any_flags(flags, BodyConfigs::OverrideFlags::LockLinearY))
            value.linear_velocity.y = real(0);
        if (math::has_any_flags(flags, BodyConfigs::OverrideFlags::LockRot))
            value.angular_velocity = real(0);
        if (!value.linear_velocity.isfinite() || !std::isfinite(value.angular_velocity))
            return stage_error(state, "advance_velocity: check failed (!value.linear_velocity.isfinite() || !std::isfinite(value.angular_velocity))");
        return state.body_pool.set_state(body, value) ||
               stage_error(state, "advance_velocity: failed to store body state");
    }
    return true;
}

bool advance_position(WorldState& state, BodyHandle body, real seconds) {
    BodyState value;
    BodyConfigs::RegisterType type;
    if (!state.body_pool.get_state(body, value) || !state.body_pool.get_register_type(body, type))
        return stage_error(state, "advance_position: check failed (!state.body_pool.get_state(body, value) || !state.body_pool.get_register_type(body, type))");
    if (type == BodyConfigs::RegisterType::ImmutableStatic)
        return true;
    integrate_body_position_seconds(value.position, value.rotation, value.linear_velocity,
                                    value.angular_velocity, seconds);
    if (!value.position.isfinite() || !std::isfinite(value.rotation))
        return stage_error(state, "advance_position: integrated body pose is not finite");
    return state.body_pool.set_state(body, value) ||
           stage_error(state, "advance_position: failed to store body state");
}

void clear_consumed(WorldState& state) {
    auto& pairs = state.pairs;
    pairs.combined_damping_ratios.clear();
    pairs.combined_contact_frequencies.clear();
    pairs.combined_max_penetration_correction_speeds.clear();
    pairs.combined_static_frictions.clear();
    pairs.combined_dynamic_frictions.clear();
    pairs.combined_restitutions.clear();
    auto& solve = state.solve;
    solve.islands.clear();
    solve.bodies.clear();
    solve.body_start_times.clear();
    solve.component_parents.clear();
    solve.component_starts.clear();
    solve.body_step_seconds.clear();
    solve.constraint_step_seconds.clear();
    solve.prepared.clear();
    solve.constraints.clear();
    solve.constraint_pair_indices.clear();
}

std::size_t component_root(WorldState::SolveData& solve, std::size_t position) {
    while (solve.component_parents[position] != position) {
        solve.component_parents[position] =
            solve.component_parents[solve.component_parents[position]];
        position = solve.component_parents[position];
    }
    return position;
}

bool prepare_components(WorldState& state, tick_t target, u32 steps) {
    auto& solve = state.solve;
    auto& pairs = state.pairs;
    const std::size_t count = solve.islands.size();
    solve.component_parents.resize(count);
    solve.component_starts.resize(count);
    if (++solve.island_position_epoch == 0) {
        std::fill(solve.island_position_stamps.begin(),
                  solve.island_position_stamps.end(), 0);
        solve.island_position_epoch = 1;
    }
    for (std::size_t i = 0; i < count; ++i) {
        const IslandId id = solve.islands[i];
        if (id.index >= solve.island_positions.size()) {
            solve.island_positions.resize(std::size_t(id.index) + 1, count);
            solve.island_position_stamps.resize(std::size_t(id.index) + 1, 0);
        }
        solve.island_positions[id.index] = i;
        solve.island_position_stamps[id.index] = solve.island_position_epoch;
        solve.component_parents[i] = i;
        solve.component_starts[i] = state.island_pool.get(id)->time;
    }
    if (!state.interval.microstep) {
        for (std::size_t row = 0; row < pairs.pairs->size(); ++row) {
            if (pairs.manifolds->counts[row] == 0)
                continue;
            const Pair& pair = (*pairs.pairs)[row];
            const auto a = state.body_pool.island_of(pair.body_a);
            const auto b = state.body_pool.island_of(pair.body_b);
            BodyConfigs::RegisterType type_a, type_b;
            if (!state.body_pool.get_register_type(pair.body_a, type_a) ||
                !state.body_pool.get_register_type(pair.body_b, type_b))
                return stage_error(state, "prepare_components: check failed (!state.body_pool.get_register_type(pair.body_a, type_a) || !state.body_pool.get_register_type(pair.body_b, type_b))");
            if (type_a == BodyConfigs::RegisterType::ImmutableStatic ||
                type_b == BodyConfigs::RegisterType::ImmutableStatic || !a || !b || *a == *b)
                continue;
            if (a->index >= solve.island_position_stamps.size() ||
                b->index >= solve.island_position_stamps.size() ||
                solve.island_position_stamps[a->index] != solve.island_position_epoch ||
                solve.island_position_stamps[b->index] != solve.island_position_epoch)
                return stage_error(state, "prepare_components: check failed (a->index >= solve.island_position_stamps.size() || b->index >= solve.island_position_stamps.size() || solve.island_position_stamps[a->index] != solve.island_position_epoch || solve.island_position_stamps[b->index] != solve.island_position_epoch)");
            const std::size_t ai = solve.island_positions[a->index];
            const std::size_t bi = solve.island_positions[b->index];
            if (ai >= count || bi >= count)
                return stage_error(state, "prepare_components: check failed (ai >= count || bi >= count)");
            const std::size_t ar = component_root(solve, ai);
            const std::size_t br = component_root(solve, bi);
            if (ar != br) {
                const std::size_t low = std::min(ar, br);
                const std::size_t high = std::max(ar, br);
                solve.component_parents[high] = low;
                solve.component_starts[low] =
                    std::max(solve.component_starts[low], solve.component_starts[high]);
            }
        }
    }
    solve.body_step_seconds.resize(solve.bodies.size());
    for (std::size_t i = 0; i < solve.bodies.size(); ++i) {
        const auto island = state.body_pool.island_of(solve.bodies[i]);
        if (!island || island->index >= solve.island_positions.size() ||
            solve.island_position_stamps[island->index] != solve.island_position_epoch)
            return stage_error(state, "prepare_components: check failed (!island || island->index >= solve.island_positions.size() || solve.island_position_stamps[island->index] != solve.island_position_epoch)");
        const std::size_t position = solve.island_positions[island->index];
        if (position >= count)
            return stage_error(state, "prepare_components: check failed (position >= count)");
        const tick_t start = state.interval.microstep
                                 ? solve.body_start_times[i]
                                 : solve.component_starts[component_root(solve, position)];
        if (start > target)
            return stage_error(state, "prepare_components: check failed (start > target)");
        const real duration = math::tu_to_seconds<real>(target - start);
        const real step = duration / static_cast<real>(steps);
        if (!std::isfinite(step))
            return stage_error(state, "prepare_components: check failed (!std::isfinite(step))");
        solve.body_step_seconds[i] = step;
    }
    return true;
}

} // namespace

bool run_solve(WorldState& state) {
    const tick_t target =
        state.interval.microstep ? state.interval.microstep_time : state.interval.target_end;
    auto& pairs = state.pairs;
    auto& solve = state.solve;
    if (++solve.selection_epoch == 0) {
        std::fill(solve.island_stamps.begin(), solve.island_stamps.end(), 0);
        for (auto& stamps : solve.body_stamps)
            std::fill(stamps.begin(), stamps.end(), 0);
        solve.selection_epoch = 1;
    }
    if (!valid_layout(pairs))
        return stage_error(state, "contact pair and manifold columns have invalid layout");
    if (!state.settings.gravity.isfinite())
        return stage_error(state, "gravity is not finite");
    if (!std::isfinite(state.settings.contact.restitution_threshold_m_s) ||
        state.settings.contact.restitution_threshold_m_s < real(0))
        return stage_error(state, "restitution threshold is invalid");

    solve.processed_islands.clear();
    solve.islands.clear();
    solve.bodies.clear();
    solve.body_start_times.clear();
    solve.component_parents.clear();
    solve.component_starts.clear();
    solve.body_step_seconds.clear();
    solve.constraint_step_seconds.clear();
    solve.constraints.clear();
    solve.prepared.clear();
    solve.constraint_pair_indices.clear();
    if (!state.interval.microstep) {
        if (!state.interval.islands)
            return stage_error(state, "run_solve: check failed (!state.interval.islands)");
        for (IslandId id : *state.interval.islands)
            if (!add_island(state, id, target))
                return stage_error(state, "run_solve: check failed (!add_island(state, id, target))");
    }
    for (std::size_t i = 0; i < pairs.pairs->size(); ++i) {
        const Pair& pair = (*pairs.pairs)[i];
        if (pair.body_a == pair.body_b || !state.body_pool.is_valid(pair.body_a) ||
            !state.body_pool.is_valid(pair.body_b) ||
            !valid_manifold(pairs.manifolds->at(i)))
            return stage_error(state, "run_solve: check failed (pair.body_a == pair.body_b || !state.body_pool.is_valid(pair.body_a) || !state.body_pool.is_valid(pair.body_b) || !valid_manifold(pairs.manifolds->at(i)))");
#ifdef SQP_DEBUG_VALIDATION
        if (!valid_parameters(pairs, i))
            return stage_error(state, "run_solve: check failed (!valid_parameters(pairs, i))");
#endif
        const auto a = state.body_pool.island_of(pair.body_a);
        const auto b = state.body_pool.island_of(pair.body_b);
        BodyConfigs::RegisterType type_a, type_b;
        if (!state.body_pool.get_register_type(pair.body_a, type_a) ||
            !state.body_pool.get_register_type(pair.body_b, type_b) ||
            (!a && type_a != BodyConfigs::RegisterType::ImmutableStatic) ||
            (!b && type_b != BodyConfigs::RegisterType::ImmutableStatic))
            return stage_error(state, "run_solve: check failed (!state.body_pool.get_register_type(pair.body_a, type_a) || !state.body_pool.get_register_type(pair.body_b, type_b) || (!a && type_a != BodyConfigs::RegisterType::ImmutableStatic) || (!b && type_b != BodyConfigs::RegisterType::ImmutableStatic))");
        if (pairs.manifolds->counts[i] > 0 &&
            ((a && state.island_pool.is_sleeping(*a) && !wake_island(state, *a)) ||
             (b && state.island_pool.is_sleeping(*b) && !wake_island(state, *b))))
            return stage_error(state, "run_solve: check failed (pairs.manifolds->counts[i] > 0 && ((a && state.island_pool.is_sleeping(*a) && !wake_island(state, *a)) || (b && state.island_pool.is_sleeping(*b) && !wake_island(state, *b))))");
        const bool enabled = pairs.manifolds->counts[i] > 0;
        if ((a && (enabled || !state.island_pool.is_sleeping(*a)) &&
             !add_island(state, *a, target)) ||
            (b && (enabled || !state.island_pool.is_sleeping(*b)) &&
             !add_island(state, *b, target)))
            return stage_error(state, "run_solve: check failed ((a && (enabled || !state.island_pool.is_sleeping(*a)) && !add_island(state, *a, target)) || (b && (enabled || !state.island_pool.is_sleeping(*b)) && !add_island(state, *b, target)))");
        if (pairs.manifolds->counts[i] > 0) {
            const tick_t start = std::max(a ? state.island_pool.get(*a)->time : tick_t(0),
                                          b ? state.island_pool.get(*b)->time : tick_t(0));
            if (!std::isfinite(math::tu_to_seconds<real>(target - start)))
                return stage_error(state, "run_solve: check failed (!std::isfinite(math::tu_to_seconds<real>(target - start)))");
        }
    }
    if (!collect_bodies(state, target))
        return stage_error(state, "run_solve: check failed (!collect_bodies(state, target))");
    for (IslandId id : solve.islands)
        solve.processed_islands.push_back({id, state.island_pool.get(id)->time, target});
    for (std::size_t i = 0; i < solve.bodies.size(); ++i)
        if (!valid_body(state, solve.bodies[i], solve.body_start_times[i], target))
            return stage_error(state, "run_solve: check failed (!valid_body(state, solve.bodies[i], solve.body_start_times[i], target))");

    const u32 steps =
        state.interval.microstep ? 1u : std::clamp(state.settings.contact.substeps, 1u, 64u);
    if (!prepare_components(state, target, steps))
        return stage_error(state, "run_solve: check failed (!prepare_components(state, target, steps))");

    const std::size_t count = pairs.pairs->size();
    state.post_solve.reports.clear();
    solve.constraints.reserve(count);
    solve.constraint_pair_indices.reserve(count);
    solve.prepared.reserve(count);

    if (!state.interval.microstep) {
        for (std::size_t i = 0; i < solve.bodies.size(); ++i) {
            const auto island = state.body_pool.island_of(solve.bodies[i]);
            const std::size_t position = solve.island_positions[island->index];
            const tick_t start = solve.component_starts[component_root(solve, position)];
            if (start == solve.body_start_times[i])
                continue;
            const real seconds = math::tu_to_seconds<real>(start - solve.body_start_times[i]);
            if (!advance_velocity(state, solve.bodies[i], seconds) ||
                !advance_position(state, solve.bodies[i], seconds))
                return stage_error(state, "run_solve: check failed (!advance_velocity(state, solve.bodies[i], seconds) || !advance_position(state, solve.bodies[i], seconds))");
        }
    }
    for (std::size_t i = 0; i < solve.bodies.size(); ++i)
        if (!advance_velocity(state, solve.bodies[i], solve.body_step_seconds[i]))
            return stage_error(state, "run_solve: check failed (!advance_velocity(state, solve.bodies[i], solve.body_step_seconds[i]))");
    // A TOI impulse acts at the impact pose, after motion with the incoming velocity.
    if (state.interval.microstep)
        for (std::size_t i = 0; i < solve.bodies.size(); ++i)
            if (!advance_position(state, solve.bodies[i], solve.body_step_seconds[i]))
                return stage_error(state, "run_solve: check failed (!advance_position(state, solve.bodies[i], solve.body_step_seconds[i]))");

    for (std::size_t i = 0; i < count; ++i) {
        const ContactManifold manifold = pairs.manifolds->at(i);
        if (manifold.count == 0)
            continue;
        const Pair& pair = (*pairs.pairs)[i];
        const auto a = state.body_pool.island_of(pair.body_a);
        const auto b = state.body_pool.island_of(pair.body_b);
        tick_t start = std::max(a ? state.island_pool.get(*a)->time : tick_t(0),
                                b ? state.island_pool.get(*b)->time : tick_t(0));
        if (!state.interval.microstep) {
            BodyConfigs::RegisterType type_a;
            if (!state.body_pool.get_register_type(pair.body_a, type_a))
                return stage_error(state, "run_solve: check failed (!state.body_pool.get_register_type(pair.body_a, type_a))");
            const auto island = type_a == BodyConfigs::RegisterType::ImmutableStatic ? b : a;
            if (island) {
                const std::size_t position = solve.island_positions[island->index];
                start = solve.component_starts[component_root(solve, position)];
            }
        }
        if (start == target)
            continue;
        const ContactSolverParameters parameters{
            pairs.combined_restitutions[i],
            pairs.combined_dynamic_frictions[i],
            pairs.combined_static_frictions[i],
            pairs.combined_damping_ratios[i],
            pairs.combined_contact_frequencies[i],
            pairs.combined_max_penetration_correction_speeds[i]};
        const real h = math::tu_to_seconds<real>(target - start) / static_cast<real>(steps);
        if (!(h > real(0)) || !std::isfinite(h))
            return stage_error(state, "run_solve: check failed (!(h > real(0)) || !std::isfinite(h))");
        const bool local = state.interval.microstep &&
                           pairs.toi_local_contacts.size() == pairs.pairs->size() &&
                           pairs.toi_local_contacts[i].valid;
        auto constraint = local
                              ? prepare_toi_local_contact_constraint(
                                    state.body_pool, pair.body_a, pair.body_b,
                                    pairs.toi_local_contacts[i].manifold,
                                    pairs.toi_local_contacts[i].origin, parameters, h, false,
                                    state.settings.contact.restitution_threshold_m_s)
                              : prepare_contact_constraint(
                                    state.body_pool, pair.body_a, pair.body_b, manifold, parameters,
                                    h, false, state.settings.contact.restitution_threshold_m_s);
        if (constraint.count != manifold.count)
            return stage_error(state, "run_solve: check failed (constraint.count != manifold.count)");
        restore_contact_impulses(solve, pair,
                                 local ? pairs.toi_local_contacts[i].manifold : manifold,
                                 start, h, constraint);
        solve.constraints.push_back(std::move(constraint));
        solve.constraint_pair_indices.push_back(i);
        solve.constraint_step_seconds.push_back(h);
    }
    solve.prepared.resize(solve.constraints.size());
    for (std::size_t i = 0; i < solve.constraints.size(); ++i)
        if (!prepare_contact_runtime(state.body_pool, solve.constraints[i], solve.prepared[i]))
            return stage_error(state, "run_solve: check failed (!prepare_contact_runtime(state.body_pool, solve.constraints[i], solve.prepared[i]))");

    const u32 iterations = std::max(1u, state.settings.contact.velocity_iterations);
    // Allow up to penetration_slop_m of overlap before applying penetration-correction bias.
    constexpr real penetration_slop_m = real(0);
    for (u32 substep = 0; substep < steps; ++substep) {
        if (substep != 0)
            for (std::size_t i = 0; i < solve.bodies.size(); ++i)
                if (!advance_velocity(state, solve.bodies[i], solve.body_step_seconds[i]))
                    return stage_error(state, "run_solve: check failed (!advance_velocity(state, solve.bodies[i], solve.body_step_seconds[i]))");
        for (std::size_t i = 0; i < solve.constraints.size(); ++i) {
            auto& constraint = solve.constraints[i];
            const real h = solve.constraint_step_seconds[i];
            const std::size_t row = solve.constraint_pair_indices[i];
            constraint.inverse_h = real(1) / h;
            const real hz = pairs.combined_contact_frequencies[row];
            const bool fixed_a =
                std::max(constraint.inverse_mass_a.x, constraint.inverse_mass_a.y) == real(0);
            const bool fixed_b =
                std::max(constraint.inverse_mass_b.x, constraint.inverse_mass_b.y) == real(0);
            constraint.softness = make_contact_softness((fixed_a || fixed_b) ? real(2) * hz : hz,
                                                        pairs.combined_damping_ratios[row], h);
        }
        for (auto& prepared : solve.prepared)
            warm_start_contact(prepared);
        for (u32 pass = 0; pass < iterations; ++pass)
            for (auto& prepared : solve.prepared)
                solve_contact_constraint(prepared, penetration_slop_m, true);
        if (!state.interval.microstep)
            for (std::size_t i = 0; i < solve.bodies.size(); ++i)
                if (!advance_position(state, solve.bodies[i], solve.body_step_seconds[i]))
                    return stage_error(state, "run_solve: check failed (!advance_position(state, solve.bodies[i], solve.body_step_seconds[i]))");
        for (auto& prepared : solve.prepared)
            refresh_contact_constraint_geometry(prepared);
        for (u32 pass = 0; pass < iterations; ++pass)
            for (auto& prepared : solve.prepared)
                solve_contact_constraint(prepared, real(0), false);
    }
    for (IslandId id : solve.islands)
        if (!state.island_pool.set_time(id, target))
            return stage_error(state, "run_solve: check failed (!state.island_pool.set_time(id, target))");
    for (std::size_t i = 0; i < solve.prepared.size(); ++i) {
        auto& prepared = solve.prepared[i];
        apply_contact_restitution(prepared);
        const std::size_t row = solve.constraint_pair_indices[i];
        const Pair& pair = (*pairs.pairs)[row];
        const auto& constraint = *prepared.constraint;
        if (state.interval.microstep &&
            std::any_of(constraint.points, constraint.points + constraint.count,
                        [](const ContactConstraintPoint& point) {
                            return point.normal_impulse > real(0) || point.tangent_impulse != real(0);
                        })) {
            for (const auto body : {pair.body_a, pair.body_b})
                if (const auto island = state.body_pool.island_of(body))
                    solve.processed_islands[solve.island_positions[island->index]]
                        .received_toi_impulse = true;
        }
        const bool local = state.interval.microstep &&
                           pairs.toi_local_contacts.size() == pairs.pairs->size() &&
                           pairs.toi_local_contacts[row].valid;
        save_contact_impulses(solve, pair,
                              local ? pairs.toi_local_contacts[row].manifold
                                    : pairs.manifolds->at(row),
                              *prepared.constraint, target,
                              solve.constraint_step_seconds[i]);
        BodyConfigs::EventFlags flags_a, flags_b;
        if (!state.body_pool.get_event_flags(pair.body_a, flags_a) ||
            !state.body_pool.get_event_flags(pair.body_b, flags_b))
            return stage_error(state, "run_solve: check failed (!state.body_pool.get_event_flags(pair.body_a, flags_a) || !state.body_pool.get_event_flags(pair.body_b, flags_b))");
        constexpr auto flag = BodyConfigs::EventFlags::AfterSolve;
        if ((flags_a & flag) != BodyConfigs::EventFlags(0) ||
            (flags_b & flag) != BodyConfigs::EventFlags(0))
            state.post_solve.reports.push_back(
                {row, report_contact_constraint(*prepared.constraint)});
    }
    for (std::size_t i = 0; i < solve.bodies.size(); ++i)
        if (solve.body_start_times[i] < target)
            bump_motion_version(state, solve.bodies[i]);
    if (!run_update_contact_edges(state))
        return stage_error(state, "run_solve: check failed (!run_update_contact_edges(state))");
    if (state.interval.microstep) {
        for (const auto& processed : solve.processed_islands)
            if (processed.received_toi_impulse &&
                !state.island_pool.increment_toi_truncations(processed.id))
                return stage_error(state, "run_solve: failed to charge TOI impulse to island");
    }
    prune_contact_cache(state);
    clear_consumed(state);
    return true;
}

bool run_advance_islands(WorldState& state, std::span<const IslandId> islands, tick_t target) {
    auto& solve = state.solve;
    if (++solve.selection_epoch == 0) {
        std::fill(solve.island_stamps.begin(), solve.island_stamps.end(), 0);
        for (auto& stamps : solve.body_stamps)
            std::fill(stamps.begin(), stamps.end(), 0);
        solve.selection_epoch = 1;
    }
    solve.islands.clear();
    solve.bodies.clear();
    solve.body_start_times.clear();
    if (!state.settings.gravity.isfinite())
        return stage_error(state, "run_advance_islands: check failed (!state.settings.gravity.isfinite())");
    for (IslandId id : islands)
        if (state.island_pool.frozen_until(id) > target) {
            continue;
        } else if (state.island_pool.is_sleeping(id)) {
            if (!state.island_pool.set_time(id, target))
                return stage_error(state, "run_advance_islands: check failed (!state.island_pool.set_time(id, target))");
        } else if (!add_island(state, id, target))
            return stage_error(state, "run_advance_islands: check failed (!add_island(state, id, target))");
    if (!collect_bodies(state, target))
        return stage_error(state, "run_advance_islands: check failed (!collect_bodies(state, target))");
    for (std::size_t i = 0; i < solve.bodies.size(); ++i)
        if (!valid_body(state, solve.bodies[i], solve.body_start_times[i], target))
            return stage_error(state, "run_advance_islands: check failed (!valid_body(state, solve.bodies[i], solve.body_start_times[i], target))");
    for (std::size_t i = 0; i < solve.bodies.size(); ++i)
        if (!advance_velocity(state, solve.bodies[i],
                              math::tu_to_seconds<real>(target - solve.body_start_times[i])))
            return stage_error(state, "run_advance_islands: check failed (!advance_velocity(state, solve.bodies[i], math::tu_to_seconds<real>(target - solve.body_start_times[i])))");
    for (std::size_t i = 0; i < solve.bodies.size(); ++i)
        if (!advance_position(state, solve.bodies[i],
                              math::tu_to_seconds<real>(target - solve.body_start_times[i])))
            return stage_error(state, "run_advance_islands: check failed (!advance_position(state, solve.bodies[i], math::tu_to_seconds<real>(target - solve.body_start_times[i])))");
    for (IslandId id : solve.islands)
        if (!state.island_pool.set_time(id, target))
            return stage_error(state, "run_advance_islands: check failed (!state.island_pool.set_time(id, target))");
    for (std::size_t i = 0; i < solve.bodies.size(); ++i)
        if (solve.body_start_times[i] < target)
            bump_motion_version(state, solve.bodies[i]);
    solve.islands.clear();
    solve.bodies.clear();
    solve.body_start_times.clear();
    prune_contact_cache(state);
    return true;
}

} // namespace sqp
