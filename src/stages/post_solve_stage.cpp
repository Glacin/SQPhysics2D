// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "world/world_stages.h"

#include <cmath>
#include <limits>

namespace sqp {
namespace {

bool valid_manifold_layout(const ContactManifolds& m, std::size_t count) {
    if (count > std::numeric_limits<std::size_t>::max() / 2)
        return false;
    const std::size_t points = count * 2;
    return m.counts.size() == count && m.normals_x.size() == count && m.normals_y.size() == count &&
           m.keys.size() == points && m.points_x.size() == points && m.points_y.size() == points &&
           m.depths.size() == points;
}

bool valid_manifold(const ContactManifolds& m, std::size_t row) {
    const u8 count = m.counts[row];
    if (count > 2)
        return false;
#ifdef SQP_DEBUG_VALIDATION
    if (!std::isfinite(m.normals_x[row]) || !std::isfinite(m.normals_y[row]))
        return false;
    for (std::size_t point = 0; point < count; ++point) {
        const std::size_t slot = row * 2 + point;
        if (!std::isfinite(m.points_x[slot]) || !std::isfinite(m.points_y[slot]) ||
            !std::isfinite(m.depths[slot]))
            return false;
    }
#endif
    return true;
}

#ifdef SQP_DEBUG_VALIDATION
bool valid_report(const ContactSolveReport& report) {
    if (!std::isfinite(report.normal_impulse) || !std::isfinite(report.tangent_impulse) ||
        !std::isfinite(report.total_impulse))
        return false;
    for (std::size_t point = 0; point < 2; ++point) {
        if (!std::isfinite(report.points[point].normal_impulse) ||
            !std::isfinite(report.points[point].tangent_impulse))
            return false;
    }
    return true;
}
#endif

void clear_consumed(WorldState& state) {
    auto& pairs = state.pairs;
    pairs.original_manifolds.clear();
    pairs.original_manifolds_valid = false;
    state.post_solve.reports.clear();
    state.post_solve.events.clear();
}

} // namespace

bool run_post_solve(WorldState& state) {
    const auto& pairs = state.pairs;
    const auto& reports = state.post_solve.reports;
    if (!pairs.pairs || !pairs.manifolds)
        return stage_error(state, "contact pair or manifold storage is missing");
    if (!pairs.original_manifolds_valid)
        return stage_error(state, "original contact manifolds are unavailable");

    const std::size_t count = pairs.pairs->size();
    if (!valid_manifold_layout(*pairs.manifolds, count) ||
        !valid_manifold_layout(pairs.original_manifolds, count) ||
        count > std::numeric_limits<std::size_t>::max() / 2)
        return stage_error(state, "run_post_solve: check failed (!valid_manifold_layout(*pairs.manifolds, count) || !valid_manifold_layout(pairs.original_manifolds, count) || count > std::numeric_limits<std::size_t>::max() / 2)");

    for (std::size_t row = 0; row < count; ++row) {
        if (!valid_manifold(*pairs.manifolds, row) ||
            !valid_manifold(pairs.original_manifolds, row))
            return stage_error(state, "run_post_solve: check failed (!valid_manifold(*pairs.manifolds, row) || !valid_manifold(pairs.original_manifolds, row))");
#ifdef SQP_DEBUG_VALIDATION
        const Pair& pair = (*pairs.pairs)[row];
        BodyConfigs::EventFlags flags_a, flags_b;
        BodyHandle owner_a, owner_b;
        if (pair.body_a == pair.body_b || pair.shape_a == pair.shape_b ||
            !state.body_pool.get_event_flags(pair.body_a, flags_a) ||
            !state.body_pool.get_event_flags(pair.body_b, flags_b) ||
            !state.body_pool.get_shape_owner(pair.shape_a, owner_a) ||
            !state.body_pool.get_shape_owner(pair.shape_b, owner_b) || owner_a != pair.body_a ||
            owner_b != pair.body_b)
            return stage_error(state, "run_post_solve: contact pair has invalid body, shape, event flags, or owner");
#endif
    }
#ifdef SQP_DEBUG_VALIDATION
    std::size_t previous = 0;
    bool first = true;
#endif
    for (const auto& entry : reports) {
        const std::size_t row = entry.pair_index;
        if (row >= count || pairs.manifolds->counts[row] == 0)
            return stage_error(state, "run_post_solve: check failed (row >= count || pairs.manifolds->counts[row] == 0)");
#ifdef SQP_DEBUG_VALIDATION
        if ((!first && row <= previous) || !valid_report(entry.contact))
            return stage_error(state, "run_post_solve: check failed ((!first && row <= previous) || !valid_report(entry.contact))");
        first = false;
        previous = row;
#endif
    }

    auto& events = state.post_solve.events;
    events.clear();
    const auto callback = state.callbacks.post_solve;
    const auto context = state.callbacks.post_data;
    if (callback) {
        events.reserve(reports.size());
        for (const auto& entry : reports) {
            const std::size_t row = entry.pair_index;
            events.push_back(
                {(*pairs.pairs)[row], pairs.original_manifolds.at(row), entry.contact});
        }
    }
    for (const PostSolveEvent& contact : events) {
        PhysicsEvent event{PhysicsEventTypes::PostSolve, contact};
        callback(event, context);
    }
    clear_consumed(state);
    return true;
}

} // namespace sqp
