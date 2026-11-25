#include "solver/integration.h"
#include "world/world_stages.h"

#include <algorithm>
#include <cmath>
#include <limits>
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
#include <chrono>
#endif

namespace sqp {
namespace {

bool valid_bounds(const Aabb& bounds) {
    return bounds.is_valid() && std::isfinite(bounds.min.x) && std::isfinite(bounds.min.y) &&
           std::isfinite(bounds.max.x) && std::isfinite(bounds.max.y);
}

bool prepare_sweep(const WorldState& state, BodyHandle body, const Island* island,
                   BodyState& start_state, BodyState& end_state, u32& flags) {
    BodyConfigs::RegisterType registration;
    if (!state.body_pool.get_state(body, start_state) ||
        !state.body_pool.get_register_type(body, registration))
        return false;
    const bool static_body = registration == BodyConfigs::RegisterType::ImmutableStatic;
    const auto id = state.body_pool.island_of(body);
    const bool frozen = id && state.island_pool.frozen_until(*id) > state.schedule.time;
    if ((!static_body && !island) ||
        (island && island->time > state.interval.target_end && !frozen))
        return false;
    end_state = start_state;
    if (!static_body && !frozen)
        predict_body_state(start_state.position, start_state.rotation, start_state.linear_velocity,
                           start_state.angular_velocity, {}, real(0),
                           state.interval.target_end - island->time, end_state.position,
                           end_state.rotation);
    flags = static_body ? BroadPhase::PROXY_STATIC
                        : registration == BodyConfigs::RegisterType::Kinematic
                              ? BroadPhase::PROXY_KINEMATIC : BroadPhase::PROXY_DYNAMIC;
    return true;
}

bool compute_bounds(const WorldState& state, ShapeHandle shape, const BodyState& start_state,
                    const BodyState& end_state, Aabb& start, Aabb& end, Aabb& swept) {
    start = state.shape_pool.world_aabb(shape, start_state);
    end = state.shape_pool.world_aabb(shape, end_state);
    swept = Aabb::merge(start, end);
    if (math::abs(end_state.rotation - start_state.rotation) > config::EPSILON) {
        const real radius = state.shape_pool.sweep_radius(shape);
        const Aabb rotational_bounds{
            {std::min(start_state.position.x, end_state.position.x) - radius,
             std::min(start_state.position.y, end_state.position.y) - radius},
            {std::max(start_state.position.x, end_state.position.x) + radius,
             std::max(start_state.position.y, end_state.position.y) + radius}};
        swept = Aabb::merge(swept, rotational_bounds);
    }
    return valid_bounds(start) && valid_bounds(end) && valid_bounds(swept);
}

bool candidate_allowed(const WorldState::BroadPhaseData& broad, BroadPhase::Pair pair) {
    if (pair.a >= broad.bodies.size() || pair.b >= broad.bodies.size() || pair.a == pair.b)
        return false;
    const auto& a = broad.filter_configurations[pair.a];
    const auto& b = broad.filter_configurations[pair.b];
    return !math::has_any_flags(a.configs.override_flags,
                                BodyConfigs::OverrideFlags::NoCollision) &&
           !math::has_any_flags(b.configs.override_flags,
                                BodyConfigs::OverrideFlags::NoCollision) &&
           !(math::has_any_flags(a.configs.override_flags,
                                 BodyConfigs::OverrideFlags::IsLiteBullet) &&
             math::has_any_flags(b.configs.override_flags,
                                 BodyConfigs::OverrideFlags::IsLiteBullet)) &&
           should_collide(a.collision_filter, b.collision_filter);
}

} // namespace

bool run_broad_phase(WorldState& state) {
    const auto& interval = state.interval;
    auto& broad = state.broad_phase;
    const auto& active_islands = broad.reprepare_after_toi
                                     ? broad.toi_changed_islands : *interval.islands;
    broad.tree.set_proxy_extension_ratio(state.settings.broad_phase.proxy_extension_ratio);
    const auto update_proxy = [&](i32 handle, const Aabb& bounds, u32 flags) {
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
        const auto begin = std::chrono::steady_clock::now();
        bool reinserted = false;
        const bool updated = broad.tree.update_proxy(handle, bounds, flags, &reinserted);
        state.stats.pending_update_nanoseconds += static_cast<u64>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - begin).count());
        state.stats.pending_proxy_reinsertions += reinserted ? 1 : 0;
        return updated;
#else
        return broad.tree.update_proxy(handle, bounds, flags);
#endif
    };

    if (broad.shapes.size() != broad.bodies.size() ||
        broad.start_aabbs.size() != broad.bodies.size() ||
        broad.end_aabbs.size() != broad.bodies.size() ||
        broad.proxy_handles.size() != broad.bodies.size())
        return stage_error(state, "proxy body, shape, AABB, and handle columns have inconsistent sizes");

    u32 epoch = broad.selection_epoch;
    if (broad.reprepare_after_toi || !active_islands.empty()) {
        if (++broad.selection_epoch == 0) {
            std::fill(broad.selected.begin(), broad.selected.end(), 0);
            broad.selection_epoch = 1;
        }
        epoch = broad.selection_epoch;
    }

    auto& interval_indices = broad.interval_proxy_indices;
    interval_indices.clear();
    auto& selected = broad.selected;
    selected.resize(broad.bodies.size(), 0);

    if (!broad.proxy_index_ready && !broad.bodies.empty()) {
        // Direct stage fixtures populate proxy columns without the command-side index.
        for (IslandId id : active_islands) {
            const Island* island = state.island_pool.get(id);
            if (!island || (island->time > interval.target_end &&
                            state.island_pool.frozen_until(id) <= state.schedule.time))
                return stage_error(state, "run_broad_phase: check failed (!island || island->time > interval.target_end)");
        }
        broad.filter_configurations.resize(broad.bodies.size());
        broad.proxy_suppressed.resize(broad.bodies.size(), 0);
        for (std::size_t i = 0; i < broad.bodies.size(); ++i) {
            const BodyHandle body = broad.bodies[i];
            if (!state.body_pool.is_valid(body) ||
                !state.body_pool.get_configuration(body, broad.filter_configurations[i]))
                return stage_error(state, "run_broad_phase: check failed (!state.body_pool.is_valid(body) || !state.body_pool.get_configuration(body, broad.filter_configurations[i]))");
            auto& handle = broad.proxy_handles[i];
            auto& suppressed = broad.proxy_suppressed[i];
            const bool disabled = math::has_any_flags(
                broad.filter_configurations[i].configs.override_flags,
                BodyConfigs::OverrideFlags::NoCollision);
            if (disabled) {
                if (handle >= 0) {
                    if (!broad.tree.remove_proxy(handle)) return stage_error(state, "run_broad_phase: check failed (!broad.tree.remove_proxy(handle))");
                    handle = -1;
                }
                suppressed = 1;
                continue;
            }
            const bool restore = handle == -1 && suppressed != 0;
            if (!restore && !broad.tree.proxy_matches(handle, static_cast<u32>(i)))
                return stage_error(state, "run_broad_phase: check failed (!restore && !broad.tree.proxy_matches(handle, static_cast<u32>(i)))");
            suppressed = 0;
            const auto id = state.body_pool.island_of(body);
            const bool in_interval = id &&
                std::find(active_islands.begin(), active_islands.end(), *id) != active_islands.end();
            if (!in_interval && !restore) continue;
            const Island* island = id ? state.island_pool.get(*id) : nullptr;
            Aabb start, end, swept;
            BodyState start_state, end_state;
            u32 flags;
            if (!prepare_sweep(state, body, island, start_state, end_state, flags) ||
                !compute_bounds(state, broad.shapes[i], start_state, end_state, start, end,
                                swept))
                return stage_error(state, "run_broad_phase: check failed (!prepare_sweep(state, body, island, start_state, end_state, flags) || !compute_bounds(state, broad.shapes[i], start_state, end_state, start, end, swept))");
            if (restore) {
                handle = broad.tree.add_proxy(static_cast<u32>(i), swept, flags);
                if (handle < 0) return stage_error(state, "run_broad_phase: check failed (handle < 0)");
            } else if (!update_proxy(handle, swept, flags)) {
                return stage_error(state, "run_broad_phase: check failed (!update_proxy(handle, swept, flags))");
            }
            broad.start_aabbs[i] = start;
            broad.end_aabbs[i] = end;
            if (in_interval) {
                selected[i] = epoch;
                interval_indices.push_back(static_cast<u32>(i));
            }
        }
    } else for (IslandId id : active_islands) {
        const Island* island = state.island_pool.get(id);
        if (!island || (island->time > interval.target_end &&
                        state.island_pool.frozen_until(id) <= state.schedule.time))
            return stage_error(state, "run_broad_phase: check failed (!island || island->time > interval.target_end)");
        const auto visit = [&](BodyHandle body) {
            if (broad.bodies.empty()) return true;
            const auto dense = state.body_pool.dense_index(body);
            if (!dense) return stage_error(state, "run_broad_phase: check failed (!dense)");
            const auto& by_kind = broad.proxies_by_dense[static_cast<std::size_t>(body.kind())];
            if (*dense >= by_kind.size()) return stage_error(state, "run_broad_phase: check failed (*dense >= by_kind.size())");
            BodyState start_state, end_state;
            u32 flags = 0;
            bool sweep_prepared = false;
            for (u32 i : by_kind[*dense]) {
                if (i >= broad.bodies.size() || broad.bodies[i] != body ||
                    !state.shape_pool.is_valid(broad.shapes[i]))
                    return stage_error(state, "run_broad_phase: check failed (i >= broad.bodies.size() || broad.bodies[i] != body || !state.shape_pool.is_valid(broad.shapes[i]))");
                if (broad.proxy_suppressed[i])
                    continue;
                if (!sweep_prepared) {
                    if (!prepare_sweep(state, body, island, start_state, end_state, flags))
                        return stage_error(state, "run_broad_phase: check failed (!prepare_sweep(state, body, island, start_state, end_state, flags))");
                    sweep_prepared = true;
                }
                Aabb start, end, swept;
                if (!compute_bounds(state, broad.shapes[i], start_state, end_state, start, end,
                                    swept) ||
                    !update_proxy(broad.proxy_handles[i], swept, flags))
                    return stage_error(state, "run_broad_phase: check failed (!compute_bounds(state, broad.shapes[i], start_state, end_state, start, end, swept) || !update_proxy(broad.proxy_handles[i], swept, flags))");
                broad.start_aabbs[i] = start;
                broad.end_aabbs[i] = end;
                selected[i] = epoch;
                interval_indices.push_back(i);
            }
            return true;
        };
        if (const auto* small = std::get_if<Island::SmallIsland>(&island->bodies)) {
            for (BodyHandle body : small->data) {
                if (!is_valid(body)) break;
                if (!visit(body)) return stage_error(state, "run_broad_phase: check failed (!visit(body))");
            }
        } else {
            for (BodyHandle body : std::get<std::vector<BodyHandle>>(island->bodies))
                if (!visit(body)) return stage_error(state, "run_broad_phase: check failed (!visit(body))");
        }
    }

    if (broad.reprepare_after_toi) {
        broad.candidates = broad.ordinary_candidates;
        std::erase_if(broad.candidates, [&](BroadPhase::Pair pair) {
            return (pair.a < selected.size() && selected[pair.a] == epoch) ||
                   (pair.b < selected.size() && selected[pair.b] == epoch);
        });
    } else {
        broad.ordinary_candidates.clear();
    }
    std::erase_if(broad.candidates,
                  [&](BroadPhase::Pair pair) { return !candidate_allowed(broad, pair); });
    const std::size_t retained_count = broad.candidates.size();
    // Marks suppress the second query of a pair whose two owners are both in this interval.
    for (u32 index : interval_indices) {
        const std::size_t previous = broad.candidates.size();
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
        BroadPhase::QueryStats query_stats;
        const auto begin = std::chrono::steady_clock::now();
        const bool queried = broad.tree.query_proxy_pairs(
            broad.proxy_handles[index], broad.candidates, &broad.traversal_stack,
            &query_stats, selected, epoch);
        state.stats.pending_query_nanoseconds += static_cast<u64>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - begin).count());
        state.stats.pending_nodes_visited += query_stats.nodes_visited;
        ++state.stats.pending_query_calls;
        state.stats.freeze_retry_queries += broad.freeze_retry;
#else
        const bool queried = broad.tree.query_proxy_pairs(
            broad.proxy_handles[index], broad.candidates, &broad.traversal_stack,
            nullptr, selected, epoch);
#endif
        if (!queried)
            return stage_error(state, "run_broad_phase: check failed (!queried)");
        auto first_new = broad.candidates.begin() + previous;
        broad.candidates.erase(
            std::remove_if(first_new, broad.candidates.end(),
                           [&](BroadPhase::Pair pair) { return !candidate_allowed(broad, pair); }),
            broad.candidates.end());
    }
    // Fat bounds change tree topology. Keep newly found pairs in a stable solver order.
    std::sort(broad.candidates.begin() + retained_count, broad.candidates.end());
    if (retained_count != 0)
        std::inplace_merge(broad.candidates.begin(),
                           broad.candidates.begin() + retained_count,
                           broad.candidates.end());
    broad.ordinary_candidates = broad.candidates;
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
    const auto rebuild_begin = std::chrono::steady_clock::now();
#endif
    broad.tree.rebuild_changed();
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
    state.stats.pending_rebuild_nanoseconds += static_cast<u64>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - rebuild_begin).count());
#endif
    return true;
}

} // namespace sqp
