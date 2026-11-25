// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "world/world_stages.h"

#include <algorithm>
#include <cmath>
#include <type_traits>

namespace sqp {
namespace {

bool rebuild_proxies(WorldState& state) {
    auto& broad = state.broad_phase;
    broad.tree.clear();
    broad.tree.set_proxy_extension_ratio(state.settings.broad_phase.proxy_extension_ratio);
    broad.bodies.clear();
    broad.proxy_index_ready = false;
    for (auto& by_kind : broad.proxies_by_dense) by_kind.clear();
    broad.shapes.clear();
    broad.start_aabbs.clear();
    broad.end_aabbs.clear();
    broad.proxy_handles.clear();
    broad.proxy_suppressed.clear();
    broad.filter_configurations.clear();
    broad.candidate_query_scratch.clear();
    broad.interval_proxy_indices.clear();
    broad.toi_changed_islands.clear();
    broad.reprepare_after_toi = false;
    broad.freeze_retry = false;
    state.narrow_phase.prepared_time.reset();
    state.narrow_phase.pair_results.clear();
    state.narrow_phase.ordinary_contact_keys.clear();
    broad.selected.clear();
    broad.candidates.clear();
    broad.ordinary_candidates.clear();
    state.narrow_phase.candidate_cache.clear();
    state.narrow_phase.toi_pairs.clear();
    state.narrow_phase.toi_times.clear();
    state.narrow_phase.toi_normals.clear();
    state.narrow_phase.toi_impacts.clear();

    // Predicted contacts were computed before these body/shape changes.
    state.pairs.pairs->clear();
    state.pairs.manifolds->clear();
    state.pairs.toi_local_contacts.clear();
    state.pairs.original_manifolds.clear();
    state.pairs.original_manifolds_valid = false;
    state.pairs.combined_damping_ratios.clear();
    state.pairs.combined_contact_frequencies.clear();
    state.pairs.combined_max_penetration_correction_speeds.clear();
    state.pairs.combined_static_frictions.clear();
    state.pairs.combined_dynamic_frictions.clear();
    state.pairs.combined_restitutions.clear();
    state.post_solve.reports.clear();
    state.post_solve.events.clear();

    std::vector<BodyHandle> bodies;
    std::vector<ShapeHandle> shapes;
    state.body_pool.get_handles(bodies);
    for (BodyConfigs::Kind kind :
         {BodyConfigs::Kind::Rigid, BodyConfigs::Kind::LiteCircle, BodyConfigs::Kind::AabbBox})
        broad.proxies_by_dense[static_cast<std::size_t>(kind)].resize(
            state.body_pool.configs_column(kind).size());
    for (BodyHandle body : bodies) {
        BodyState pose;
        BodyConfiguration configuration;
        if (!state.body_pool.get_state(body, pose) ||
            !state.body_pool.get_configuration(body, configuration) ||
            !state.body_pool.get_attached_shapes(body, shapes))
            return stage_error(state, "rebuild_proxies: check failed (!state.body_pool.get_state(body, pose) || !state.body_pool.get_configuration(body, configuration) || !state.body_pool.get_attached_shapes(body, shapes))");
        for (ShapeHandle shape : shapes) {
            if (!state.shape_pool.is_valid(shape) ||
                broad.bodies.size() >= std::numeric_limits<u32>::max())
                return stage_error(state, "rebuild_proxies: check failed (!state.shape_pool.is_valid(shape) || broad.bodies.size() >= std::numeric_limits<u32>::max())");
            const Aabb bounds = state.shape_pool.world_aabb(shape, pose);
            if (!bounds.is_valid() || !bounds.min.isfinite() || !bounds.max.isfinite())
                return stage_error(state, "rebuild_proxies: check failed (!bounds.is_valid() || !bounds.min.isfinite() || !bounds.max.isfinite())");
            const u32 index = static_cast<u32>(broad.bodies.size());
            const bool suppressed = math::has_any_flags(configuration.configs.override_flags,
                                                        BodyConfigs::OverrideFlags::NoCollision);
            const u32 flags =
                configuration.configs.register_type == BodyConfigs::RegisterType::ImmutableStatic
                    ? BroadPhase::PROXY_STATIC
                    : configuration.configs.register_type == BodyConfigs::RegisterType::Kinematic
                          ? BroadPhase::PROXY_KINEMATIC : BroadPhase::PROXY_DYNAMIC;
            const i32 proxy = suppressed ? -1 : broad.tree.add_proxy(index, bounds, flags);
            if (!suppressed && proxy < 0)
                return stage_error(state, "rebuild_proxies: check failed (!suppressed && proxy < 0)");
            broad.bodies.push_back(body);
            const auto dense = state.body_pool.dense_index(body);
            if (!dense) return stage_error(state, "rebuild_proxies: check failed (!dense)");
            broad.proxies_by_dense[static_cast<std::size_t>(body.kind())][*dense].push_back(index);
            broad.shapes.push_back(shape);
            broad.start_aabbs.push_back(bounds);
            broad.end_aabbs.push_back(bounds);
            broad.proxy_handles.push_back(proxy);
            broad.proxy_suppressed.push_back(suppressed ? 1 : 0);
            broad.filter_configurations.push_back(configuration);
        }
    }
    broad.proxy_index_ready = true;
    return true;
}

bool update_proxies(WorldState& state, const std::vector<BodyHandle>& changed) {
    auto& broad = state.broad_phase;
    broad.toi_changed_islands.clear();
    broad.reprepare_after_toi = false;
    broad.freeze_retry = false;
    state.narrow_phase.prepared_time.reset();
    state.narrow_phase.pair_results.clear();
    state.narrow_phase.ordinary_contact_keys.clear();
    broad.candidates.clear();
    broad.ordinary_candidates.clear();
    broad.interval_proxy_indices.clear();
    state.narrow_phase.toi_pairs.clear();
    state.narrow_phase.toi_times.clear();
    state.narrow_phase.toi_normals.clear();
    state.narrow_phase.toi_impacts.clear();
    state.pairs.pairs->clear();
    state.pairs.manifolds->clear();
    state.pairs.toi_local_contacts.clear();
    state.pairs.original_manifolds.clear();
    state.pairs.original_manifolds_valid = false;
    state.pairs.combined_damping_ratios.clear();
    state.pairs.combined_contact_frequencies.clear();
    state.pairs.combined_max_penetration_correction_speeds.clear();
    state.pairs.combined_static_frictions.clear();
    state.pairs.combined_dynamic_frictions.clear();
    state.pairs.combined_restitutions.clear();
    state.post_solve.reports.clear();
    state.post_solve.events.clear();
    for (BodyHandle body : changed) {
        if (!state.body_pool.is_valid(body)) continue;
        const auto dense = state.body_pool.dense_index(body);
        if (!dense) return stage_error(state, "update_proxies: check failed (!dense)");
        const auto& by_kind = broad.proxies_by_dense[static_cast<std::size_t>(body.kind())];
        if (*dense >= by_kind.size()) return stage_error(state, "update_proxies: check failed (*dense >= by_kind.size())");
        BodyState pose;
        BodyConfiguration config;
        if (!state.body_pool.get_state(body, pose) ||
            !state.body_pool.get_configuration(body, config))
            return stage_error(state, "update_proxies: check failed (!state.body_pool.get_state(body, pose) || !state.body_pool.get_configuration(body, config))");
        const bool suppressed = math::has_any_flags(
            config.configs.override_flags, BodyConfigs::OverrideFlags::NoCollision);
        const u32 flags = config.configs.register_type ==
                                  BodyConfigs::RegisterType::ImmutableStatic
                              ? BroadPhase::PROXY_STATIC
                              : config.configs.register_type ==
                                        BodyConfigs::RegisterType::Kinematic
                                    ? BroadPhase::PROXY_KINEMATIC : BroadPhase::PROXY_DYNAMIC;
        for (u32 index : by_kind[*dense]) {
            if (index >= broad.bodies.size() || broad.bodies[index] != body ||
                !state.shape_pool.is_valid(broad.shapes[index]))
                return stage_error(state, "update_proxies: check failed (index >= broad.bodies.size() || broad.bodies[index] != body || !state.shape_pool.is_valid(broad.shapes[index]))");
            const Aabb bounds = state.shape_pool.world_aabb(broad.shapes[index], pose);
            if (!bounds.is_valid() || !bounds.min.isfinite() || !bounds.max.isfinite())
                return stage_error(state, "update_proxies: check failed (!bounds.is_valid() || !bounds.min.isfinite() || !bounds.max.isfinite())");
            auto& handle = broad.proxy_handles[index];
            if (suppressed) {
                if (handle >= 0 && !broad.tree.remove_proxy(handle)) return stage_error(state, "update_proxies: check failed (handle >= 0 && !broad.tree.remove_proxy(handle))");
                handle = -1;
            } else if (handle < 0) {
                handle = broad.tree.add_proxy(index, bounds, flags);
                if (handle < 0) return stage_error(state, "update_proxies: check failed (handle < 0)");
            } else if (!broad.tree.update_proxy(handle, bounds, flags)) {
                return stage_error(state, "update_proxies: check failed (!broad.tree.update_proxy(handle, bounds, flags))");
            }
            broad.proxy_suppressed[index] = suppressed ? 1 : 0;
            broad.filter_configurations[index] = config;
            broad.start_aabbs[index] = bounds;
            broad.end_aabbs[index] = bounds;
        }
    }
    return true;
}

bool register_body(WorldState& state, BodyHandle body, tick_t now) {
    if (!is_valid(body))
        return stage_error(state, "register_body: check failed (!is_valid(body))");
    BodyConfigs::RegisterType registration;
    if (!state.body_pool.get_register_type(body, registration))
        return stage_error(state, "register_body: check failed (!state.body_pool.get_register_type(body, registration))");
    if (registration == BodyConfigs::RegisterType::ImmutableStatic)
        return true;
    const IslandId island = state.island_pool.create();
    if (!island.valid()) {
        state.body_pool.remove(body);
        return stage_error(state, "register_body: check failed (!island.valid())");
    }
    if (!state.island_pool.set_time(island, now) ||
        !state.island_pool.add_body(island, body, state.body_pool)) {
        state.island_pool.destroy(island, state.body_pool);
        state.body_pool.remove(body);
        return stage_error(state, "register_body: check failed (!state.island_pool.set_time(island, now) || !state.island_pool.add_body(island, body, state.body_pool))");
    }
    return reschedule_island(state, island);
}

CommandResult apply(WorldState& state, PhysicsCommand& command, tick_t now, bool& dirty) {
    CommandResult result;
    auto& bodies = state.body_pool;
    auto& shapes = state.shape_pool;
    bool ok = false;
    BodyHandle affected_before;
    if (const auto* removal = std::get_if<RemoveShape>(&command))
        bodies.get_shape_owner(removal->shape, affected_before);
    std::visit(
        [&](auto& c) {
            using T = std::decay_t<decltype(c)>;
            if constexpr (std::is_same_v<T, SetVelocity>) {
                ok = c.linear.isfinite() && std::isfinite(c.angular) &&
                     bodies.set_velocity(c.body, c.linear, c.angular);
                dirty |= ok;
            } else if constexpr (std::is_same_v<T, SetPose>) {
                ok = c.position.isfinite() && std::isfinite(c.rotation) &&
                     bodies.set_pose(c.body, c.position, c.rotation);
                dirty |= ok;
            } else if constexpr (std::is_same_v<T, SetBodyConfig>) {
                std::visit(
                    [&](auto value) {
                        using V = decltype(value);
                        if constexpr (std::is_same_v<V, BodyConfigs::CcdOverride>)
                            ok = bodies.set_ccd_override(c.body, value);
                        else if constexpr (std::is_same_v<V, BodyConfigs::EventFlags>)
                            ok = bodies.set_event_flags(c.body, value);
                        else
                            ok = bodies.set_override_flags(c.body, value);
                    },
                    c.value);
                dirty |= ok;
            } else if constexpr (std::is_same_v<T, RemoveBody>) {
                state.connectivity.full_rebuild = true;
                const auto island = bodies.island_of(c.body);
                ok = bodies.remove(c.body, state.island_pool);
                if (ok && island) {
                    const Island* remaining = state.island_pool.get(*island);
                    bool empty = false;
                    if (remaining) {
                        if (const auto* small =
                                std::get_if<Island::SmallIsland>(&remaining->bodies))
                            empty = !is_valid(small->data[0]);
                        else
                            empty = std::get<std::vector<BodyHandle>>(remaining->bodies).empty();
                    }
                    if (empty) {
                        ok = unschedule_island(state, *island) &&
                             state.island_pool.destroy(*island);
                        if (ok) {
                            std::erase(*state.interval.islands, *island);
                            std::erase(*state.interval.delayed_islands, *island);
                        }
                    } else {
                        ok = wake_island(state, *island);
                    }
                }
                dirty |= ok;
            } else if constexpr (std::is_same_v<T, ApplyImpulse>) {
                BodyState velocity;
                BodyProperties properties;
                BodyConfiguration config;
                if (c.impulse.isfinite() && std::isfinite(c.angular_impulse) &&
                    bodies.get_state(c.body, velocity) &&
                    bodies.get_properties(c.body, properties) &&
                    bodies.get_configuration(c.body, config) &&
                    config.configs.register_type == BodyConfigs::RegisterType::Dynamic &&
                    properties.inv_mass > 0) {
                    velocity.linear_velocity += c.impulse * properties.inv_mass;
                    velocity.angular_velocity += c.angular_impulse * properties.inv_inertia;
                    ok = bodies.set_velocity(c.body, velocity.linear_velocity,
                                             velocity.angular_velocity);
                }
                dirty |= ok;
            } else if constexpr (std::is_same_v<T, SetConstantForce>) {
                ok = bodies.set_constant_force(c.body, c.force, c.torque);
                dirty |= ok;
            } else if constexpr (std::is_same_v<T, ClearForces>) {
                ok = bodies.clear_forces(c.body);
                dirty |= ok;
            } else if constexpr (std::is_same_v<T, SetCollisionFilter>) {
                ok = bodies.set_collision_filter(c.body, c.filter);
                dirty |= ok;
            } else if constexpr (std::is_same_v<T, AttachShape>) {
                ok = bodies.attach_shape(c.body, c.shape, shapes);
                dirty |= ok;
            } else if constexpr (std::is_same_v<T, DetachShape>) {
                ok = bodies.detach_shape(c.body, c.shape, shapes);
                dirty |= ok;
            } else if constexpr (std::is_same_v<T, CreateBody>) {
                std::visit(
                    [&](auto& def) {
                        using D = std::remove_cvref_t<decltype(*def)>;
                        if (!def)
                            return;
                        if constexpr (std::is_same_v<D, RigidBodyDef>)
                            result.body = bodies.add_rigid(*def);
                        else if constexpr (std::is_same_v<D, LiteCircleBodyDef>)
                            result.body = bodies.add_lite_circle(*def, c.shape, shapes);
                        else
                            result.body = bodies.add_aabb_box(*def, c.shape, shapes);
                    },
                    c.definition);
                ok = register_body(state, result.body, now);
                if (!ok)
                    result.body = {};
                dirty |= ok;
            } else if constexpr (std::is_same_v<T, CreateShape>) {
                if (is_valid(c.body) &&
                    (!bodies.is_valid(c.body) || bodies.kind(c.body) != BodyConfigs::Kind::Rigid))
                    return;
                std::visit(
                    [&](const auto& def) {
                        using D = std::decay_t<decltype(def)>;
                        if constexpr (std::is_same_v<D, CircleShapeDef>)
                            result.shape = shapes.add_circle(def);
                        else
                            result.shape = shapes.add_polygon(def);
                    },
                    c.definition);
                ok = is_valid(result.shape);
                if (ok && is_valid(c.body)) {
                    ok = bodies.attach_shape(c.body, result.shape, shapes);
                    if (!ok) {
                        BodyHandle owner;
                        if (bodies.get_shape_owner(result.shape, owner) && owner == c.body)
                            bodies.detach_shape(c.body, result.shape, shapes);
                        shapes.remove(result.shape);
                        result.shape = {};
                    }
                }
                dirty |= ok;
            } else if constexpr (std::is_same_v<T, RemoveShape>) {
                BodyHandle owner;
                const bool attached = bodies.get_shape_owner(c.shape, owner);
                ok = shapes.is_valid(c.shape) &&
                     (!attached || bodies.detach_shape(owner, c.shape, shapes)) &&
                     shapes.remove(c.shape);
                dirty |= ok;
            } else if constexpr (std::is_same_v<T, SetSettings>) {
                state.connectivity.full_rebuild = true;
                state.settings = c.settings;
                state.settings.sanitize();
                ok = true;
                dirty = true;
                auto sleeping = state.schedule.period_islands[0];
                for (IslandId id : sleeping)
                    ok &= wake_island(state, id);
            } else if constexpr (std::is_same_v<T, SetCallback>) {
                if (c.type == PhysicsEventTypes::PreSolve) {
                    state.callbacks.pre_solve = c.callback;
                    state.callbacks.pre_data = c.callback ? c.user_data : nullptr;
                    ok = true;
                } else if (c.type == PhysicsEventTypes::PostSolve) {
                    state.callbacks.post_solve = c.callback;
                    state.callbacks.post_data = c.callback ? c.user_data : nullptr;
                    ok = true;
                }
            }
        },
        command);
    if (ok && dirty) {
        state.connectivity.dirty = true;
        BodyHandle touched = affected_before;
        std::visit(
            [&](const auto& c) {
                using T = std::decay_t<decltype(c)>;
                if constexpr (requires { c.body; })
                    touched = c.body;
                else if constexpr (std::is_same_v<T, RemoveShape>)
                    bodies.get_shape_owner(c.shape, touched);
            },
            command);
        if (bodies.is_valid(touched)) {
            if (const auto dense = bodies.dense_index(touched)) {
                auto versions = bodies.motion_version_column(touched.kind());
                if (*dense < versions.size())
                    ++versions[*dense];
            }
            if (const auto island = bodies.island_of(touched)) {
                state.connectivity.dirty_islands.push_back(*island);
                ok &= wake_island(state, *island);
            }
        }
    }
    result.status = ok ? CommandStatus::Applied : CommandStatus::Failed;
    return result;
}

} // namespace

bool run_process_command(WorldState& state) {
    const auto& interval = state.interval;
    const tick_t now = interval.microstep ? interval.microstep_time : interval.target_end;
    auto& commands = state.commands;
    bool dirty = false;
    bool structural = false;
    std::vector<BodyHandle> changed_bodies;
    while (!commands.heap.empty() && commands.heap.front().time <= now) {
        std::pop_heap(commands.heap.begin(), commands.heap.end(), WorldState::CommandData::Later{});
        auto pending = std::move(commands.heap.back());
        commands.heap.pop_back();
        CommandResult result;
        result.status = CommandStatus::Expired;
        if (pending.time == now) {
            auto& affected = state.schedule.scratch_islands;
            affected.clear();
            if (std::holds_alternative<SetSettings>(pending.command)) {
                state.island_pool.get_ids(affected);
            } else {
                BodyHandle body;
                std::visit(
                    [&](const auto& c) {
                        using T = std::decay_t<decltype(c)>;
                        if constexpr (requires { c.body; })
                            body = c.body;
                        else if constexpr (std::is_same_v<T, RemoveShape>)
                            state.body_pool.get_shape_owner(c.shape, body);
                    },
                    pending.command);
                if (const auto island = state.body_pool.island_of(body))
                    affected.push_back(*island);
            }
            std::erase_if(affected, [&](IslandId id) {
                u8 bucket;
                std::size_t position;
                return state.island_pool.frozen_until(id) > now ||
                       (state.island_pool.get_schedule_location(id, bucket, position) && bucket == 0);
            });
            if (!run_advance_islands(state, affected, now))
                return stage_error(state, "run_process_command: check failed (!run_advance_islands(state, affected, now))");
            result = apply(state, pending.command, now, dirty);
            // A failed command is a command result, not a stage failure.
            if (result.status == CommandStatus::Failed)
                state.stage_error_detail.clear();
            if (result.status == CommandStatus::Applied) {
                std::visit([&](const auto& c) {
                    using T = std::decay_t<decltype(c)>;
                    if constexpr (std::is_same_v<T, CreateBody> ||
                                  std::is_same_v<T, RemoveBody> ||
                                  std::is_same_v<T, CreateShape> ||
                                  std::is_same_v<T, RemoveShape> ||
                                  std::is_same_v<T, AttachShape> ||
                                  std::is_same_v<T, DetachShape> ||
                                  std::is_same_v<T, SetSettings>) {
                        structural = true;
                    } else if constexpr (requires { c.body; }) {
                        changed_bodies.push_back(c.body);
                    }
                }, pending.command);
            }
        }
        commands.results.push_back({pending.id, result});
    }
    if (!dirty)
        return true;
    const bool updated = structural ? rebuild_proxies(state) : update_proxies(state, changed_bodies);
    return updated || stage_error(state, "run_process_command: failed to synchronize broad-phase proxies");
}

} // namespace sqp
