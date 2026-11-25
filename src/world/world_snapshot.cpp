// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "sqphysics2d/detail/core_world.h"

#include <stdexcept>
#include <type_traits>
#include <utility>

namespace sqp {
namespace {

PhysicsCommand copy_command(const PhysicsCommand& command) {
    return std::visit([](const auto& value) -> PhysicsCommand {
        using T = std::decay_t<decltype(value)>;
        if constexpr (std::is_same_v<T, CreateBody>) {
            CreateBody copy;
            copy.shape = value.shape;
            copy.definition = std::visit([](const auto& definition) -> decltype(copy.definition) {
                using Definition = std::remove_cv_t<std::remove_reference_t<decltype(*definition)>>;
                return definition ? std::make_unique<Definition>(*definition)
                                  : std::unique_ptr<Definition>{};
            }, value.definition);
            return copy;
        } else {
            return value;
        }
    }, command);
}

} // namespace

WorldState::WorldState(const WorldState& source)
    : stage_error_detail(source.stage_error_detail),
      last_error_message(source.last_error_message), stats(source.stats),
      settings(source.settings), body_pool(source.body_pool), shape_pool(source.shape_pool),
      island_pool(source.island_pool), contact_graph(source.contact_graph),
      connectivity(source.connectivity), broad_phase(source.broad_phase),
      narrow_phase(source.narrow_phase), process_toi(source.process_toi),
      schedule(source.schedule), solve(source.solve), post_solve(source.post_solve),
      callbacks(source.callbacks) {
    joint_pool.copy_from(source.joint_pool, body_pool);

    commands.next_id = source.commands.next_id;
    commands.results = source.commands.results;
    commands.heap.reserve(source.commands.heap.size());
    for (const auto& pending : source.commands.heap)
        commands.heap.push_back({pending.time, pending.id, copy_command(pending.command)});

    *pairs.pairs = *source.pairs.pairs;
    *pairs.manifolds = *source.pairs.manifolds;
    pairs.toi_local_contacts = source.pairs.toi_local_contacts;
    pairs.original_manifolds = source.pairs.original_manifolds;
    pairs.original_manifolds_valid = source.pairs.original_manifolds_valid;
    *pairs.delayed_pairs = *source.pairs.delayed_pairs;
    *pairs.delayed_manifolds = *source.pairs.delayed_manifolds;
    pairs.delayed_toi_local_contacts = source.pairs.delayed_toi_local_contacts;
    pairs.delayed_original_manifolds = source.pairs.delayed_original_manifolds;
    pairs.delayed_original_manifolds_valid = source.pairs.delayed_original_manifolds_valid;
    pairs.combined_damping_ratios = source.pairs.combined_damping_ratios;
    pairs.combined_contact_frequencies = source.pairs.combined_contact_frequencies;
    pairs.combined_max_penetration_correction_speeds =
        source.pairs.combined_max_penetration_correction_speeds;
    pairs.combined_static_frictions = source.pairs.combined_static_frictions;
    pairs.combined_dynamic_frictions = source.pairs.combined_dynamic_frictions;
    pairs.combined_restitutions = source.pairs.combined_restitutions;

    *interval.islands = *source.interval.islands;
    *interval.delayed_islands = *source.interval.delayed_islands;
    interval.microstep = source.interval.microstep;
    interval.tick = source.interval.tick;
    interval.microstep_time = source.interval.microstep_time;
    interval.target_end = source.interval.target_end;

    process_toi.toi_keys.clear();
    for (const auto& event : process_toi.toi_heap)
        process_toi.toi_keys.insert({event.time, event.pair});
    solve.contact_cache_index.clear();
    for (std::size_t i = 0; i < solve.contact_cache.size(); ++i)
        solve.contact_cache_index.emplace(detail::contact_pair_key(solve.contact_cache[i]), i);

    for (auto& entry : narrow_phase.candidate_cache)
        if (entry.view.pool)
            entry.view.pool = &shape_pool;
    for (std::size_t i = 0; i < solve.prepared.size(); ++i) {
        if (solve.prepared[i].constraint && i < solve.constraints.size()) {
            if (!prepare_contact_runtime(body_pool, solve.constraints[i], solve.prepared[i]))
                throw std::runtime_error("WorldState copy: invalid prepared contact binding");
        } else {
            solve.prepared[i] = {};
        }
    }
}

WorldSnapshot::WorldSnapshot() = default;
WorldSnapshot::~WorldSnapshot() = default;
WorldSnapshot::WorldSnapshot(const WorldSnapshot& other)
    : state_(other.state_ ? std::make_unique<WorldState>(*other.state_) : nullptr) {}
WorldSnapshot& WorldSnapshot::operator=(const WorldSnapshot& other) {
    if (this != &other) {
        WorldSnapshot copy(other);
        state_.swap(copy.state_);
    }
    return *this;
}
WorldSnapshot::WorldSnapshot(WorldSnapshot&&) noexcept = default;
WorldSnapshot& WorldSnapshot::operator=(WorldSnapshot&&) noexcept = default;

} // namespace sqp
