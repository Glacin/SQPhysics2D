#pragma once

#include "sqphysics2d/detail/core_world.h"

#include <bit>
#include <type_traits>

namespace sqp::detail {

// FNV-1a over explicitly encoded scalar fields. Never feed object representations:
// padding, capacity, and pointer values are not part of a simulation state.
class StateHash {
public:
    void byte(u8 value) noexcept {
        value_ ^= value;
        value_ *= 1099511628211ull;
    }
    void integer(u64 value) noexcept {
        for (unsigned i = 0; i < 8; ++i) {
            byte(static_cast<u8>(value));
            value >>= 8;
        }
    }
    template <class T> void add(T value) noexcept
        requires(std::is_integral_v<T> || std::is_enum_v<T>)
    {
        if constexpr (std::is_enum_v<T>)
            integer(static_cast<u64>(static_cast<std::underlying_type_t<T>>(value)));
        else
            integer(static_cast<u64>(value));
    }
    void add(float value) noexcept { integer(std::bit_cast<u32>(value)); }
    void add(double value) noexcept { integer(std::bit_cast<u64>(value)); }
    void add(Vec2 value) noexcept { add(value.x); add(value.y); }
    void add(Aabb value) noexcept { add(value.min); add(value.max); }
    void add(BodyHandle value) noexcept { add(value.kind()); add(value.index()); add(value.generation); }
    void add(ShapeHandle value) noexcept { add(value.kind); add(value.index); add(value.generation); }
    void add(IslandId value) noexcept { add(value.index); add(value.generation); }
    void add(JointHandle value) noexcept { add(value.index); add(value.generation); }
    void add(const BodyConfigs& value) noexcept {
        add(value.override_flags); add(value.runtime_flags); add(value.event_flags);
        add(value.ccd_override); add(value.register_type); add(value.kind);
    }
    void add(CollisionFilter value) noexcept { add(value.category_bits); add(value.collision_mask); }
    void add(const PhysicsMaterial& value) noexcept {
        add(value.contact_frequency_hz); add(value.max_penetration_correction_speed);
        add(value.contact_damping_ratio); add(value.base_restitution);
        add(value.restitution_impact_threshold); add(value.restitution_impact_exponent);
        add(value.dynamic_friction); add(value.static_friction);
        add(value.restitution_priority); add(value.friction_priority);
        add(value.contact_damping_priority); add(value.contact_response_priority);
        add(value.restitution_flags);
    }
    void add(const BodyState& value) noexcept {
        add(value.position); add(value.rotation); add(value.linear_velocity);
        add(value.angular_velocity); add(value.flags);
    }
    void add(const BodyProperties& value) noexcept {
        add(value.mass); add(value.inv_mass); add(value.inertia); add(value.inv_inertia);
        add(value.local_center); add(value.bound_radius); add(value.density);
        add(value.radius); add(value.half_extents);
    }
    void add(const BodyConfiguration& value) noexcept {
        add(value.drag_multiplier); add(value.configs); add(value.collision_filter);
    }
    void add(const ContactFeatureRef& value) noexcept {
        add(value.type); add(value.index); add(value.auxiliary);
    }
    void add(const ContactPointKey& value) noexcept { add(value.feature_a); add(value.feature_b); }
    void add(const Pair& value) noexcept {
        add(value.body_a); add(value.body_b); add(value.shape_a); add(value.shape_b);
        add(value.motion_version_a); add(value.motion_version_b);
    }
    void add(BroadPhase::Pair value) noexcept { add(value.a); add(value.b); }
    void add(const ConstraintGraph::ContactEdge& value) noexcept {
        add(value.body_a); add(value.body_b); add(value.shape_a); add(value.shape_b);
    }
    void add(const PhysicsSettings& value) noexcept {
        add(value.interval.max_linear_displacement_m);
        add(value.interval.max_angular_displacement_rad);
        add(value.interval.min_period_tu); add(value.interval.max_truncations_per_tick);
        add(value.contact.substeps); add(value.contact.velocity_iterations);
        add(value.contact.toi_velocity_iterations); add(value.contact.restitution_threshold_m_s);
        add(value.contact.enable_dynamic_restitution);
        add(value.ccd.enable_circle); add(value.ccd.enable_polygon);
        add(value.ccd.linear_start_threshold_m);
        add(value.broad_phase.proxy_extension_ratio);
        add(value.gravity); add(value.linear_velocity_tick_multiplier);
        add(value.angular_velocity_tick_multiplier); add(value.allow_sleep);
        add(value.sleep_linear_velocity_threshold_m_s);
        add(value.sleep_angular_velocity_threshold_rad_s);
        add(value.sleep_time_threshold_tu);
    }
    template <class T> void sequence(const T& values) {
        add(static_cast<u64>(values.size()));
        for (const auto& item : values)
            add(item);
    }
    u64 finish() const noexcept { return value_; }
private:
    u64 value_ = 14695981039346656037ull;
};

u64 hash_world_state(const WorldState& state);

} // namespace sqp::detail
