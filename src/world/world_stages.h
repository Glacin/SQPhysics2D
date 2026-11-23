// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#pragma once

#include "sqphysics2d/detail/core_world.h"
#include <span>

namespace sqp {

inline bool stage_error(WorldState& state, std::string_view detail) {
    if (state.stage_error_detail.empty())
        state.stage_error_detail = detail;
    return false;
}

// Consumes command-heap entries at the currently selected interval time.
bool run_process_command(WorldState& state);

// Synchronizes collision-disabled proxies, updates interval proxies, and appends filtered
// candidates. Returns false when the interval or a proxy is invalid.
bool run_broad_phase(WorldState& state);

// Classifies all broad-phase candidates, appends outputs, then consumes the candidates on success.
bool run_narrow_phase(WorldState& state);
tick_t ceil_toi_time(tick_t start, tick_t end, real fraction);

// Consumes narrow-phase TOI records and selects the earliest valid in-bound time.
// selected_time is true only when this invocation found such a time.
bool run_process_toi(WorldState& state, tick_t bound, bool& selected_time);

// Removes ordinary pairs whose body motion versions have changed, keeping the matching
// contact manifolds in the same order. Returns false for inconsistent manifold columns.
bool run_update_pairs(WorldState& state);

// Validates active contacts and writes callback-adjusted manifolds and material results atomically.
bool run_pre_solve(WorldState& state);

// Advances active islands and solves the contact constraints prepared by PreSolve.
bool run_solve(WorldState& state);
// Solve-owned no-contact integration used before time-stamped command mutations.
bool run_advance_islands(WorldState& state, std::span<const IslandId> islands, tick_t target);
bool freeze_island_until_next_tick(WorldState& state, IslandId id, tick_t now);

// Collect islands, change islands' interval, sleep on / off.
bool run_collect_tasks(WorldState& state);
bool reschedule_island(WorldState& state, IslandId id, tick_t elapsed = 0);
bool wake_island(WorldState& state, IslandId id);
bool park_island(WorldState& state, IslandId id);
bool unschedule_island(WorldState& state, IslandId id);

// Choose the next time based on Toi heap, interval buckets and command queue.
bool run_advance_interval(WorldState& state);

void clear_active_interval_data(WorldState& state);

// Update the constraint graph.
bool run_update_contact_edges(WorldState& state);

// Merge / split of islands. TODO: Check this carefully and write more test cases.
bool run_rebuild_islands(WorldState& state, bool force = true);

// Dispatches solved-contact events, then clears data consumed by the solve pipeline.
bool run_post_solve(WorldState& state);

// The microstep interval flag bypasses UpdatePairs, even without contact pairs.
bool should_run_update_pairs(const WorldState& state);

} // namespace sqp
