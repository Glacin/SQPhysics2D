#include "sqphysics2d/detail/core_world.h"

#include "world_stages.h"

#include <algorithm>
#include <cstdio>
#include <limits>
#include <utility>
#include <variant>

#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
#include <chrono>
#endif

namespace sqp {
namespace {

using Stage = WorldState::StageStats::Stage;

#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
u64 island_body_count(const Island& island) {
    if (const auto* small = std::get_if<Island::SmallIsland>(&island.bodies)) {
        u64 count = 0;
        for (BodyHandle body : small->data) {
            if (!is_valid(body))
                break;
            ++count;
        }
        return count;
    }
    return std::get<std::vector<BodyHandle>>(island.bodies).size();
}

WorldState::StageStats::IntervalEntry interval_entry(const WorldState& state) {
    WorldState::StageStats::IntervalEntry entry;
    for (IslandId id : *state.interval.islands) {
        if (const Island* island = state.island_pool.get(id)) {
            ++entry.active_islands;
            entry.active_bodies += island_body_count(*island);
        }
    }
    entry.candidates = state.stats.pending_candidates;
    entry.pairs = state.pairs.pairs->size();
    entry.query_nanoseconds = state.stats.pending_query_nanoseconds;
    entry.update_nanoseconds = state.stats.pending_update_nanoseconds;
    entry.rebuild_nanoseconds = state.stats.pending_rebuild_nanoseconds;
    entry.nodes_visited = state.stats.pending_nodes_visited;
    entry.proxy_reinsertions = state.stats.pending_proxy_reinsertions;
    entry.query_calls = state.stats.pending_query_calls;
    return entry;
}
#endif

class StageTimer {
public:
    StageTimer(WorldState& state, Stage stage)
        : entry_(state.stats.stages[static_cast<std::size_t>(stage)]) {
        ++entry_.calls;
    }

    ~StageTimer() {
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
        entry_.nanoseconds += static_cast<u64>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - start_)
                .count());
#endif
    }

private:
    WorldState::StageStats::Entry& entry_;
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
    const std::chrono::steady_clock::time_point start_ = std::chrono::steady_clock::now();
#endif
};

template <typename Fn>
bool run_stage(WorldState& state, Stage stage, Fn&& fn) {
    StageTimer timer(state, stage);
    state.stage_error_detail.clear();
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
    // AdvanceInterval clears the active columns, so read them before it runs.
    WorldState::StageStats::IntervalEntry interval;
    if (stage == Stage::AdvanceInterval)
        interval = interval_entry(state);
#endif
    const bool succeeded = std::forward<Fn>(fn)();
    if (!succeeded) {
        static constexpr std::array<const char*, static_cast<std::size_t>(Stage::Count)> names = {
            "ProcessCommand", "CollectTasks", "BroadPhase", "NarrowPhase", "ProcessToi",
            "UpdatePairs", "PreSolve", "Solve", "PostSolve", "AdvanceInterval"};
        const char* detail = state.stage_error_detail.empty()
                                 ? "stage operation failed"
                                 : state.stage_error_detail.c_str();
        state.last_error_message = std::string(names[static_cast<std::size_t>(stage)]) + ": " + detail;
        std::fprintf(stderr, "%s\n", state.last_error_message.c_str());
    }
    if (succeeded) {
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
        if (stage == Stage::BroadPhase)
            state.stats.pending_candidates += state.broad_phase.candidates.size();
#endif
        if (stage == Stage::AdvanceInterval) {
            ++state.stats.intervals;
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
            state.stats.interval_entries.push_back(interval);
            state.stats.pending_candidates = 0;
            state.stats.pending_query_nanoseconds = 0;
            state.stats.pending_update_nanoseconds = 0;
            state.stats.pending_rebuild_nanoseconds = 0;
            state.stats.pending_nodes_visited = 0;
            state.stats.pending_proxy_reinsertions = 0;
            state.stats.pending_query_calls = 0;
#endif
        }
    }
    return succeeded;
}

} // namespace

CoreWorld::~CoreWorld() = default;

CoreWorld::StageId CoreWorld::process_command(VirtualWorld& world) {
    auto& core = static_cast<CoreWorld&>(world);
    return core.process_commands() ? core.stage_id(StageRole::CollectTasks) : Error;
}

bool CoreWorld::process_commands() {
    return run_stage(world_state_, Stage::ProcessCommand,
                     [&] { return run_process_command(world_state_); });
}

// These scheduling stages are implemented with the interval scheduler.
CoreWorld::StageId CoreWorld::collect_tasks(VirtualWorld& world) {
    auto& core = static_cast<CoreWorld&>(world);
    return run_stage(core.world_state_, Stage::CollectTasks,
                     [&] { return run_collect_tasks(core.world_state_); })
               ? core.stage_id(StageRole::BroadPhase)
               : Error;
}
CoreWorld::StageId CoreWorld::advance_interval(VirtualWorld& world) {
    auto& core = static_cast<CoreWorld&>(world);
    return run_stage(core.world_state_, Stage::AdvanceInterval,
                     [&] { return run_advance_interval(core.world_state_); })
               ? Finish
               : Error;
}

bool CoreWorld::step_until(tick_t target) {
    auto& state = world_state_;
    auto& schedule = state.schedule;
    if (target < schedule.time)
        return false;
    schedule.step_target = target;

    const auto prepare_ordinary = [&] {
        if (schedule.prepared || schedule.ordinary_target == std::numeric_limits<tick_t>::max())
            return true;
        clear_active_interval_data(state);
        state.interval.target_end = schedule.ordinary_target;
        if (!run_stage(state, Stage::BroadPhase, [&] { return run_broad_phase(state); }) ||
            !run_stage(state, Stage::NarrowPhase, [&] { return run_narrow_phase(state); }))
            return false;
        state.broad_phase.reprepare_after_toi = false;
        state.broad_phase.toi_changed_islands.clear();
        schedule.prepared = true;
        return true;
    };
    const auto solve_selected = [&](bool microstep) {
        state.interval.microstep = microstep;
        state.interval.tick = !microstep && state.interval.target_end % config::TU_PER_TICK == 0;
        if (!microstep &&
            !run_stage(state, Stage::UpdatePairs, [&] { return run_update_pairs(state); }))
            return false;
        return run_stage(state, Stage::PreSolve, [&] { return run_pre_solve(state); }) &&
               run_stage(state, Stage::Solve, [&] { return run_solve(state); }) &&
               run_stage(state, Stage::PostSolve, [&] { return run_post_solve(state); }) &&
               run_stage(state, Stage::AdvanceInterval,
                         [&] { return run_advance_interval(state); });
    };

    constexpr tick_t never = std::numeric_limits<tick_t>::max();
    for (;;) {
        if (schedule.ordinary_target == 0 &&
            !run_stage(state, Stage::CollectTasks, [&] { return run_collect_tasks(state); }))
            return false;
        if (!prepare_ordinary())
            return false;
        const tick_t ordinary = schedule.ordinary_target;
        const tick_t command = state.commands.heap.empty()
                                   ? never
                                   : std::max(schedule.time, state.commands.heap.front().time);
        const tick_t before_command = command == 0 ? 0 : command - 1;
        const tick_t toi_bound = std::min({target, ordinary, before_command});
        bool selected = false;
        if (command > schedule.time && toi_bound >= schedule.time &&
            !run_stage(state, Stage::ProcessToi,
                       [&] { return run_process_toi(state, toi_bound, selected); }))
            return false;
        if (!schedule.prepared)
            continue;
        if (selected) {
            if (!solve_selected(true))
                return false;
            continue;
        }
        if (command <= target && command <= ordinary) {
            const bool ordinary_at_command = command == ordinary;
            schedule.time = command;
            state.interval.target_end = command;
            if (!process_commands())
                return false;
            schedule.prepared = false;
            if (!ordinary_at_command) {
                clear_active_interval_data(state);
                schedule.ordinary_target = 0;
            } else {
                std::erase_if(*state.interval.islands,
                              [&](IslandId id) { return state.island_pool.get(id) == nullptr; });
            }
            continue;
        }
        if (ordinary <= target) {
            if (!solve_selected(false))
                return false;
            continue;
        }
        schedule.time = target;
        return true;
    }
}

CoreWorld::CoreWorld() { install_default_stages(); }

CoreWorld::CoreWorld(const WorldState& source) : world_state_(source) {
    install_default_stages();
}

void CoreWorld::install_default_stages() {
    stage_ids_[0] = install_stage(process_command, "ProcessCommand");
    stage_ids_[1] = install_stage(collect_tasks, "CollectTasks");
    stage_ids_[2] = install_stage(broad_phase, "BroadPhase");
    stage_ids_[3] = install_stage(narrow_phase, "NarrowPhase");
    stage_ids_[4] = install_stage(process_toi, "ProcessToi");
    stage_ids_[5] = install_stage(update_pairs, "UpdatePairs");
    stage_ids_[6] = install_stage(pre_solve, "PreSolve");
    stage_ids_[7] = install_stage(solve, "Solve");
    stage_ids_[8] = install_stage(post_solve, "PostSolve");
    stage_ids_[9] = install_stage(advance_interval, "AdvanceInterval");

    set_entry_stage(stage_id(StageRole::ProcessCommand));
}

CoreWorld::StageId CoreWorld::broad_phase(VirtualWorld& world) {
    auto& core = static_cast<CoreWorld&>(world);
    return run_stage(core.world_state_, Stage::BroadPhase,
                     [&] { return run_broad_phase(core.world_state_); })
               ? core.stage_id(StageRole::NarrowPhase)
               : Error;
}

CoreWorld::StageId CoreWorld::narrow_phase(VirtualWorld& world) {
    auto& core = static_cast<CoreWorld&>(world);
    if (!run_stage(core.world_state_, Stage::NarrowPhase,
                   [&] { return run_narrow_phase(core.world_state_); }))
        return Error;
    core.world_state_.schedule.prepared = true;
    return core.stage_id(StageRole::ProcessToi);
}

CoreWorld::StageId CoreWorld::process_toi(VirtualWorld& world) {
    auto& core = static_cast<CoreWorld&>(world);
    bool selected_time = false;
    if (!run_stage(core.world_state_, Stage::ProcessToi, [&] {
            return run_process_toi(core.world_state_, core.world_state_.interval.target_end,
                                   selected_time);
        }))
        return Error;
    if (selected_time)
        return core.stage_id(StageRole::ProcessCommand);
    if (!core.world_state_.schedule.prepared) {
        clear_active_interval_data(core.world_state_);
        return core.stage_id(StageRole::BroadPhase);
    }
    if (core.world_state_.interval.microstep)
        return core.stage_id(StageRole::PreSolve);
    return core.stage_id(StageRole::UpdatePairs);
}

CoreWorld::StageId CoreWorld::update_pairs(VirtualWorld& world) {
    auto& core = static_cast<CoreWorld&>(world);
    return run_stage(core.world_state_, Stage::UpdatePairs,
                     [&] { return run_update_pairs(core.world_state_); })
               ? core.stage_id(StageRole::PreSolve)
               : Error;
}

CoreWorld::StageId CoreWorld::pre_solve(VirtualWorld& world) {
    auto& core = static_cast<CoreWorld&>(world);
    return run_stage(core.world_state_, Stage::PreSolve,
                     [&] { return run_pre_solve(core.world_state_); })
               ? core.stage_id(StageRole::Solve)
               : Error;
}

CoreWorld::StageId CoreWorld::solve(VirtualWorld& world) {
    auto& core = static_cast<CoreWorld&>(world);
    return run_stage(core.world_state_, Stage::Solve,
                     [&] { return run_solve(core.world_state_); })
               ? core.stage_id(StageRole::PostSolve)
               : Error;
}

CoreWorld::StageId CoreWorld::post_solve(VirtualWorld& world) {
    auto& core = static_cast<CoreWorld&>(world);
    return run_stage(core.world_state_, Stage::PostSolve,
                     [&] { return run_post_solve(core.world_state_); })
               ? core.stage_id(StageRole::AdvanceInterval)
               : Error;
}

} // namespace sqp
