#include "state_hash.h"

#include <type_traits>
#include <variant>

namespace sqp::detail {
namespace {

void add(StateHash& h, const ContactManifolds& m) {
    h.sequence(m.keys); h.sequence(m.points_x); h.sequence(m.points_y);
    h.sequence(m.depths); h.sequence(m.counts);
    h.sequence(m.normals_x); h.sequence(m.normals_y);
}
void add(StateHash& h, const ContactManifold& m) {
    h.add(m.count); h.add(m.normal);
    for (i32 i = 0; i < m.count && i < 2; ++i) {
        h.add(m.points[i].point); h.add(m.points[i].depth); h.add(m.points[i].key);
    }
}
void add(StateHash& h, const WorldState::ToiLocalImpact& v) {
    h.add(v.valid);
    if (v.valid) {
        h.add(v.origin); h.add(v.position_a); h.add(v.position_b);
        h.add(v.rotation_a); h.add(v.rotation_b);
    }
}
void add(StateHash& h, const WorldState::ToiLocalContact& v) {
    h.add(v.valid);
    if (v.valid) { h.add(v.origin); add(h, v.manifold); }
}
template <class T> void sequence(StateHash& h, const T& values, auto&& write) {
    h.add(static_cast<u64>(values.size()));
    for (const auto& item : values) write(item);
}
void add(StateHash& h, const CircleShapeDef& v) {
    h.add(v.radius); h.add(v.local_offset); h.add(v.local_rotation);
    h.add(v.density); h.add(v.restitution_impact_scale); h.add(v.material_index);
}
void add(StateHash& h, const PolygonShapeDef& v) {
    h.sequence(v.vertices); h.add(v.vertices.overflowed());
    h.add(v.local_offset); h.add(v.local_rotation); h.add(v.density);
    h.add(v.restitution_impact_scale); h.add(v.material_index);
}
void add(StateHash& h, const RigidBodyDef& v) {
    h.add(v.call_back_ptr != nullptr); h.add(v.position); h.add(v.rotation);
    h.add(v.linear_velocity); h.add(v.angular_velocity); h.add(v.drag_multiplier);
    h.add(v.configs); h.add(v.collision_filter);
}
template <class T> void add_lite(StateHash& h, const T& v) {
    h.add(v.position); h.add(v.linear_velocity); h.add(v.drag_multiplier);
    h.add(v.configs); h.add(v.collision_filter);
}
void add(StateHash& h, const PhysicsCommand& command) {
    h.add(command.index());
    std::visit([&](const auto& v) {
        using T = std::decay_t<decltype(v)>;
        if constexpr (std::is_same_v<T, CreateBody>) {
            h.add(v.shape); h.add(v.definition.index());
            std::visit([&](const auto& ptr) {
                h.add(static_cast<bool>(ptr));
                if (ptr) {
                    using U = std::decay_t<decltype(*ptr)>;
                    if constexpr (std::is_same_v<U, RigidBodyDef>) add(h, *ptr);
                    else add_lite(h, *ptr);
                }
            }, v.definition);
        } else if constexpr (std::is_same_v<T, CreateShape>) {
            h.add(v.body); h.add(v.definition.index());
            std::visit([&](const auto& def) { add(h, def); }, v.definition);
        } else if constexpr (std::is_same_v<T, SetSettings>) {
            h.add(v.settings);
        } else if constexpr (std::is_same_v<T, SetCallback>) {
            h.add(v.type); h.add(v.callback != nullptr); h.add(v.user_data != nullptr);
        } else if constexpr (std::is_same_v<T, SetBodyConfig>) {
            h.add(v.body); h.add(v.value.index());
            std::visit([&](const auto& flag) { h.add(flag); }, v.value);
        } else {
            if constexpr (requires { v.body; }) h.add(v.body);
            if constexpr (requires { v.shape; }) h.add(v.shape);
            if constexpr (requires { v.linear; }) h.add(v.linear);
            if constexpr (requires { v.angular; }) h.add(v.angular);
            if constexpr (requires { v.position; }) h.add(v.position);
            if constexpr (requires { v.rotation; }) h.add(v.rotation);
            if constexpr (requires { v.impulse; }) h.add(v.impulse);
            if constexpr (requires { v.angular_impulse; }) h.add(v.angular_impulse);
            if constexpr (requires { v.force; }) h.add(v.force);
            if constexpr (requires { v.torque; }) h.add(v.torque);
            if constexpr (requires { v.filter; }) h.add(v.filter);
        }
    }, command);
}
} // namespace

u64 hash_world_state(const WorldState& s) {
    StateHash h;
    h.add(u64{1}); // Encoding version.
    h.add(s.settings);
    s.body_pool.append_state_hash(h);
    s.shape_pool.append_state_hash(h);
    s.island_pool.append_state_hash(h);
    s.joint_pool.append_state_hash(h);

    h.add(s.contact_graph.version());
    h.sequence(s.contact_graph.contacts());
    h.sequence(s.contact_graph.changes());
    h.add(s.connectivity.seen_graph_version);
    h.add(s.connectivity.dirty); h.add(s.connectivity.full_rebuild);
    h.sequence(s.connectivity.dirty_islands);
    h.add(static_cast<u64>(s.connectivity.last_rebuilt_body_count));

    h.add(s.commands.next_id);
    sequence(h, s.commands.heap, [&](const auto& v) {
        h.add(v.time); h.add(v.id); add(h, v.command);
    });
    sequence(h, s.commands.results, [&](const auto& v) {
        h.add(v.id); h.add(v.result.status);
        h.add(v.result.body); h.add(v.result.shape);
    });

    s.broad_phase.tree.append_state_hash(h);
    h.sequence(s.broad_phase.bodies);
    for (const auto& by_kind : s.broad_phase.proxies_by_dense)
        sequence(h, by_kind, [&](const auto& row) { h.sequence(row); });
    h.add(s.broad_phase.proxy_index_ready);
    h.sequence(s.broad_phase.shapes);
    h.sequence(s.broad_phase.start_aabbs); h.sequence(s.broad_phase.end_aabbs);
    h.sequence(s.broad_phase.proxy_handles); h.sequence(s.broad_phase.proxy_suppressed);
    h.sequence(s.broad_phase.interval_proxy_indices);
    h.sequence(s.broad_phase.selected); h.add(s.broad_phase.selection_epoch);
    h.sequence(s.broad_phase.toi_changed_islands); h.add(s.broad_phase.reprepare_after_toi);
    h.add(s.broad_phase.freeze_retry);
    h.sequence(s.broad_phase.candidates); h.sequence(s.broad_phase.ordinary_candidates);

    sequence(h, s.process_toi.toi_heap, [&](const auto& v) {
        h.add(v.time); h.add(v.pair); h.add(v.normal);
        h.add(v.sequence); add(h, v.impact);
    });
    h.add(s.process_toi.next_sequence);
    h.add(s.narrow_phase.prepared_time.has_value());
    if (s.narrow_phase.prepared_time) {
        h.add(*s.narrow_phase.prepared_time);
        h.add(s.narrow_phase.prepared_end);
    }
    h.sequence(s.narrow_phase.ordinary_contact_keys);
    h.sequence(s.narrow_phase.toi_pairs);
    h.sequence(s.narrow_phase.toi_times); h.sequence(s.narrow_phase.toi_normals);
    sequence(h, s.narrow_phase.toi_impacts, [&](const auto& v) { add(h, v); });

    h.sequence(*s.pairs.pairs); add(h, *s.pairs.manifolds);
    sequence(h, s.pairs.toi_local_contacts, [&](const auto& v) { add(h, v); });
    h.add(s.pairs.original_manifolds_valid);
    if (s.pairs.original_manifolds_valid) add(h, s.pairs.original_manifolds);
    h.sequence(*s.pairs.delayed_pairs); add(h, *s.pairs.delayed_manifolds);
    sequence(h, s.pairs.delayed_toi_local_contacts, [&](const auto& v) { add(h, v); });
    h.add(s.pairs.delayed_original_manifolds_valid);
    if (s.pairs.delayed_original_manifolds_valid) add(h, s.pairs.delayed_original_manifolds);
    h.sequence(s.pairs.combined_damping_ratios);
    h.sequence(s.pairs.combined_contact_frequencies);
    h.sequence(s.pairs.combined_max_penetration_correction_speeds);
    h.sequence(s.pairs.combined_static_frictions);
    h.sequence(s.pairs.combined_dynamic_frictions);
    h.sequence(s.pairs.combined_restitutions);

    h.sequence(*s.interval.islands); h.sequence(*s.interval.delayed_islands);
    h.add(s.interval.microstep); h.add(s.interval.tick);
    h.add(s.interval.microstep_time); h.add(s.interval.target_end);
    for (const auto& bucket : s.schedule.period_islands) h.sequence(bucket);
    sequence(h, s.schedule.entries, [&](const auto& v) {
        h.add(v.id); h.add(v.period); h.add(v.due);
        h.add(static_cast<u64>(v.position)); h.add(v.uses_ccd);
    });
    h.add(s.schedule.time); h.add(s.schedule.step_target);
    h.add(s.schedule.ordinary_target); h.add(s.schedule.prepared);
    h.add(s.schedule.last_ordinary_end.has_value());
    if (s.schedule.last_ordinary_end)
        h.add(*s.schedule.last_ordinary_end);

    sequence(h, s.solve.contact_cache, [&](const auto& v) {
        h.add(v.body_a); h.add(v.body_b); h.add(v.shape_a); h.add(v.shape_b);
        for (int i = 0; i < 2; ++i) {
            h.add(v.keys[i]); h.add(v.normal_impulses[i]); h.add(v.tangent_impulses[i]);
        }
        h.add(v.normal); h.add(v.step_seconds); h.add(v.time); h.add(v.count);
    });
    sequence(h, s.solve.processed_islands, [&](const auto& v) {
        h.add(v.id);
        h.add(v.start);
        h.add(v.end);
        h.add(v.received_toi_impulse);
    });
    h.add(s.callbacks.pre_solve != nullptr); h.add(s.callbacks.post_solve != nullptr);
    h.add(s.callbacks.pre_data != nullptr); h.add(s.callbacks.post_data != nullptr);
    return h.finish();
}
} // namespace sqp::detail
