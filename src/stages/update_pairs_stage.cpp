// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "world/world_stages.h"

#include <limits>

namespace sqp {
namespace {

bool current_version(const BodyPool& bodies, BodyHandle body, u64 expected) {
    if (!bodies.is_valid(body))
        return false;
    const auto index = bodies.dense_index(body);
    if (!index)
        return false;
    const auto versions = bodies.motion_version_column(body.kind());
    return *index < versions.size() && versions[*index] == expected;
}

bool pair_is_current(const BodyPool& bodies, const Pair& pair) {
    return current_version(bodies, pair.body_a, pair.motion_version_a) &&
           current_version(bodies, pair.body_b, pair.motion_version_b);
}

bool valid_layout(const WorldState::PairsData& narrow) {
    const auto& m = *narrow.manifolds;
    const std::size_t count = narrow.pairs->size();
    if (count > std::numeric_limits<std::size_t>::max() / 2)
        return false;
    const std::size_t points = count * 2;
    return m.counts.size() == count && m.normals_x.size() == count && m.normals_y.size() == count &&
           m.keys.size() == points && m.points_x.size() == points && m.points_y.size() == points &&
           m.depths.size() == points;
}

} // namespace

bool should_run_update_pairs(const WorldState& state) {
    return !state.interval.microstep;
}

bool run_update_pairs(WorldState& state) {
    auto& narrow = state.pairs;
    if (!valid_layout(narrow))
        return stage_error(state, "contact pair and manifold columns have invalid layout");

    auto& m = *narrow.manifolds;
    std::size_t write = 0;
    for (std::size_t read = 0; read < narrow.pairs->size(); ++read) {
        if (!pair_is_current(state.body_pool, (*narrow.pairs)[read]))
            continue;
        if (write != read) {
            (*narrow.pairs)[write] = (*narrow.pairs)[read];
            m.counts[write] = m.counts[read];
            m.normals_x[write] = m.normals_x[read];
            m.normals_y[write] = m.normals_y[read];
            for (std::size_t slot = 0; slot < 2; ++slot) {
                const std::size_t source = read * 2 + slot;
                const std::size_t target = write * 2 + slot;
                m.keys[target] = m.keys[source];
                m.points_x[target] = m.points_x[source];
                m.points_y[target] = m.points_y[source];
                m.depths[target] = m.depths[source];
            }
        }
        ++write;
    }
    narrow.pairs->resize(write);
    m.counts.resize(write);
    m.normals_x.resize(write);
    m.normals_y.resize(write);
    m.keys.resize(write * 2);
    m.points_x.resize(write * 2);
    m.points_y.resize(write * 2);
    m.depths.resize(write * 2);
    return true;
}

} // namespace sqp
