// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "world/world_stages.h"

#include <algorithm>
#include <limits>
#include <variant>

namespace sqp {

bool freeze_island_until_next_tick(WorldState& state, IslandId id, tick_t now) {
    const Island* island = state.island_pool.get(id);
    if (!island)
        return stage_error(state, "freeze_island_until_next_tick: island does not exist");
    if (state.island_pool.frozen_until(id) > now)
        return true;
    if (state.island_pool.is_sleeping(id) && !wake_island(state, id))
        return stage_error(state, "freeze_island_until_next_tick: check failed (state.island_pool.is_sleeping(id) && !wake_island(state, id))");
    const tick_t remainder = now % config::TU_PER_TICK;
    const tick_t offset = config::TU_PER_TICK - remainder;
    if (now > std::numeric_limits<tick_t>::max() - offset)
        return stage_error(state, "next tick time would overflow");
    const tick_t until = now + offset;
    if (!state.island_pool.set_time(id, until) ||
        !state.island_pool.set_frozen_until(id, until))
        return stage_error(state, "freeze_island_until_next_tick: check failed (!state.island_pool.set_time(id, until) || !state.island_pool.set_frozen_until(id, until))");
    const auto bump = [&](BodyHandle body) {
        const auto dense = state.body_pool.dense_index(body);
        if (!dense)
            return stage_error(state, "freeze_island_until_next_tick: check failed (!dense)");
        auto versions = state.body_pool.motion_version_column(body.kind());
        if (*dense >= versions.size())
            return stage_error(state, "freeze_island_until_next_tick: check failed (*dense >= versions.size())");
        ++versions[*dense];
        return true;
    };
    if (const auto* small = std::get_if<Island::SmallIsland>(&island->bodies)) {
        for (BodyHandle body : small->data) {
            if (!is_valid(body))
                break;
            if (!bump(body))
                return stage_error(state, "freeze_island_until_next_tick: check failed (!bump(body))");
        }
    } else {
        for (BodyHandle body : std::get<std::vector<BodyHandle>>(island->bodies))
            if (!bump(body))
                return stage_error(state, "freeze_island_until_next_tick: check failed (!bump(body))");
    }
    return true;
}

void clear_active_interval_data(WorldState& state) {
    state.pairs.pairs->clear();
    state.pairs.manifolds->clear();
    state.pairs.toi_local_contacts.clear();
    state.pairs.original_manifolds.clear();
    state.pairs.original_manifolds_valid = false;
    state.broad_phase.candidates.clear();
    state.narrow_phase.toi_pairs.clear();
    state.narrow_phase.toi_times.clear();
    state.narrow_phase.toi_normals.clear();
    state.narrow_phase.toi_impacts.clear();
}

bool run_advance_interval(WorldState& state) {
    auto& schedule = state.schedule;
    auto& interval = state.interval;
    if (interval.microstep) {
        schedule.time = interval.microstep_time;
        for (const auto& processed : state.solve.processed_islands) {
            const IslandId id = processed.id;
            if (state.island_pool.get(id) &&
                state.island_pool.toi_truncations(id) >=
                    state.settings.interval.max_truncations_per_tick &&
                !freeze_island_until_next_tick(state, id, schedule.time))
                return stage_error(state, "run_advance_interval: check failed (state.island_pool.get(id) && state.island_pool.toi_truncations(id) >= state.settings.interval.max_truncations_per_tick && !freeze_island_until_next_tick(state, id, schedule.time))");
        }
        clear_active_interval_data(state);
        state.pairs.pairs.swap(state.pairs.delayed_pairs);
        state.pairs.manifolds.swap(state.pairs.delayed_manifolds);
        state.pairs.toi_local_contacts.swap(state.pairs.delayed_toi_local_contacts);
        std::swap(state.pairs.original_manifolds, state.pairs.delayed_original_manifolds);
        std::swap(state.pairs.original_manifolds_valid,
                  state.pairs.delayed_original_manifolds_valid);
        interval.islands.swap(interval.delayed_islands);
        std::erase_if(*interval.islands, [&](IslandId id) {
            const Island* island = state.island_pool.get(id);
            return island && island->time > interval.target_end &&
                   state.island_pool.frozen_until(id) > schedule.time;
        });
        interval.microstep = false;
        interval.microstep_time = 0;
        interval.tick = schedule.ordinary_target != 0 &&
                        schedule.ordinary_target != std::numeric_limits<tick_t>::max() &&
                        schedule.ordinary_target % config::TU_PER_TICK == 0;
        // Only islands advanced by the microstep need new sweeps. Candidate pairs
        // for unchanged islands and future TOI records remain available.
        auto& broad = state.broad_phase;
        broad.toi_changed_islands.clear();
        for (const auto& processed : state.solve.processed_islands)
            if (const Island* island = state.island_pool.get(processed.id);
                island && island->time <= interval.target_end)
                broad.toi_changed_islands.push_back(processed.id);
        broad.reprepare_after_toi = true;
        clear_active_interval_data(state);
        state.pairs.delayed_pairs->clear();
        state.pairs.delayed_manifolds->clear();
        state.pairs.delayed_toi_local_contacts.clear();
        state.pairs.delayed_original_manifolds.clear();
        state.pairs.delayed_original_manifolds_valid = false;
        interval.delayed_islands->clear();
        schedule.prepared = false;
        return true;
    }
    schedule.time = interval.target_end;
    if (interval.tick) {
        state.island_pool.reset_toi_truncations();
        std::vector<IslandId> ids;
        state.island_pool.get_ids(ids);
        for (IslandId id : ids) {
            const tick_t until = state.island_pool.frozen_until(id);
            if (until != 0 && until <= schedule.time) {
                if (!state.island_pool.set_frozen_until(id, 0))
                    return stage_error(state, "run_advance_interval: check failed (!state.island_pool.set_frozen_until(id, 0))");
                state.connectivity.dirty = true;
                state.connectivity.full_rebuild = true;
            }
        }
    }
    clear_active_interval_data(state);
    interval.islands->clear();
    interval.tick = false;
    schedule.prepared = false;
    schedule.ordinary_target = 0;
    return true;
}

} // namespace sqp
