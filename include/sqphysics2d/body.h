// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <limits>
#include <type_traits>
#include <variant>
#include <vector>

#include "sqphysics2d/detail/math.h"

namespace sqp {

// Configuration shared by all body kinds. Runtime flags describe the stored body state;
// override and event flags control simulation behavior.
struct BodyConfigs {
    // Flags that override physics behavior.
    enum class OverrideFlags : u16 {
        // Body never goes to sleep.
        NoSleep = 1u << 0,
        // Prevent movement along X or Y.
        LockLinearX = 1u << 2,
        LockLinearY = 1u << 3,
        // Prevent rotation.
        LockRot = 1u << 4,
        LockLinear = LockLinearX | LockLinearY,
        LockAll = LockLinear | LockRot,
        // Suppress broad-phase proxies and contact candidates for this body.
        NoCollision = 1u << 5,
        // Exclude pairs where both bodies carry this flag; it does not force CCD.
        IsLiteBullet = 1u << 6,
        // Ignore gravity's X or Y component.
        NoGravityX = 1u << 8,
        NoGravityY = 1u << 9,
        NoGravity = NoGravityX | NoGravityY,
    };

    // Opt-in event flags. The current World dispatches PreSolve and AfterSolve;
    // ContactEvents and SleepEvents are stored but have no callback dispatcher yet.
    enum class EventFlags : u16 {
        // Reserved for contact begin/end notifications.
        ContactEvents = 1u << 0,
        PreSolve = 1u << 1,
        AfterSolve = 1u << 2,
        // Reserved for sleep/wake notifications.
        SleepEvents = 1u << 3,
    };

    enum class RuntimeFlags : u8 {
        IsAlive = 1u << 1,
    };

    enum class CcdOverride : u8 {
        // Choose no cast or a shape cast from motion thresholds.
        Auto = 0,
        // Force the shape-cast CCD path.
        Shapecast = 1,
        // No CCD.
        Disabled = 2
    };

    // How a body participates in the simulation.
    enum class RegisterType : u8 {
        // Fully simulated.
        Dynamic = 0,
        // Moved by velocity, unaffected by forces/contacts.
        Kinematic = 1,
        // Never moves.
        ImmutableStatic = 2,
    };

    // Storage and supported motion class of a body.
    enum class Kind : u8 {
        None = 0,
        // General rigid body with rotation, multiple shapes.
        Rigid = 1,
        // Lightweight circle body without rotation.
        LiteCircle = 2,
        // Lightweight axis-aligned box body without rotation.
        AabbBox = 3,
    };

    OverrideFlags override_flags = OverrideFlags(0);
    RuntimeFlags runtime_flags = RuntimeFlags(0);
    EventFlags event_flags = EventFlags(0);
    CcdOverride ccd_override = CcdOverride::Auto;
    RegisterType register_type = RegisterType::Dynamic;
    Kind kind = Kind::None;
};

template <class E>
concept BodyConfigFlag =
    std::is_same_v<E, BodyConfigs::OverrideFlags> || std::is_same_v<E, BodyConfigs::EventFlags> ||
    std::is_same_v<E, BodyConfigs::RuntimeFlags>;

template <BodyConfigFlag E> constexpr E operator|(E a, E b) {
    return static_cast<E>(static_cast<std::underlying_type_t<E>>(a) |
                          static_cast<std::underlying_type_t<E>>(b));
}
template <BodyConfigFlag E> constexpr E operator&(E a, E b) {
    return static_cast<E>(static_cast<std::underlying_type_t<E>>(a) &
                          static_cast<std::underlying_type_t<E>>(b));
}
template <BodyConfigFlag E> constexpr E operator~(E value) {
    return static_cast<E>(~static_cast<std::underlying_type_t<E>>(value));
}
template <BodyConfigFlag E> constexpr E& operator|=(E& a, E b) {
    return a = a | b;
}
template <BodyConfigFlag E> constexpr E& operator&=(E& a, E b) {
    return a = a & b;
}

// Bitmask collision filtering.
struct CollisionFilter {
    // Categories this body belongs to and categories it accepts for collision.
    u64 category_bits = 1;
    u64 collision_mask = ~u64(0);
};

// Two filters collide only if each one's category is accepted by the
// other's mask.
constexpr bool should_collide(const CollisionFilter& a, const CollisionFilter& b) {
    return (a.category_bits & b.collision_mask) != 0 && (b.category_bits & a.collision_mask) != 0;
}

// Applies the NoGravity* override flags to a gravity vector and returns
// the effective gravity for that body.
Vec2 gravity_with_body_overrides(Vec2, BodyConfigs::OverrideFlags);

// Current body configuration, including drag, collision filter and flags.
struct BodyConfiguration {
    real drag_multiplier = 1;
    BodyConfigs configs;
    CollisionFilter collision_filter;
};

// Creation parameters for a rigid body with rotation and multiple attachable shapes.
// Mass and inertia are derived from the shapes attached after creation.
struct RigidBodyDef {
    RigidBodyDef() = default;
    // Convenience: position, rotation, linear velocity, angular velocity.
    RigidBodyDef(const Vec2&, real rotation, const Vec2&, real angular_velocity);
    // Stored callback pointer; the current world pipeline does not invoke it.
    void (*call_back_ptr)(void*) = nullptr;
    Vec2 position;
    // Initial angle in radians.
    real rotation = 0;
    Vec2 linear_velocity;
    // Initial angular velocity in radians per second.
    real angular_velocity = 0;

    // Per-body drag multiplier.
    real drag_multiplier = 1;
    BodyConfigs configs;
    // Per-body continuous collision mode and mutual category/mask filter.
    CollisionFilter collision_filter;
};

// Creation parameters for a translation-only circle body bound to one circle shape.
struct LiteCircleBodyDef {
    LiteCircleBodyDef() = default;
    // Convenience: position, linear velocity, optional drag multiplier.
    LiteCircleBodyDef(const Vec2&, const Vec2&, real drag = 1);
    Vec2 position, linear_velocity;

    real drag_multiplier = 1;
    // Per-body continuous collision mode and mutual category/mask filter.
    BodyConfigs configs;
    CollisionFilter collision_filter;
};

// Creation parameters for a translation-only axis-aligned box bound to one polygon shape.
struct AabbBoxBodyDef {
    AabbBoxBodyDef() = default;
    // Convenience: position, linear velocity, optional drag multiplier.
    AabbBoxBodyDef(const Vec2&, const Vec2&, real drag = 1);
    Vec2 position, linear_velocity;

    real drag_multiplier = 1;
    // Per-body continuous collision mode and mutual category/mask filter.
    BodyConfigs configs;
    CollisionFilter collision_filter;
};

// Value representation of a body's pose, velocity and runtime flags.
// BodyPool stores these values in columns; lite bodies ignore rotation and angular velocity.
struct BodyState {
    Vec2 position;
    // Angle in radians.
    real rotation = 0;
    Vec2 linear_velocity;
    // Radians per second.
    real angular_velocity = 0;
    BodyConfigs::RuntimeFlags flags = BodyConfigs::RuntimeFlags(0);
    bool operator==(const BodyState&) const = default;

    // Checks if a runtime flag is set.
    bool has_flag(BodyConfigs::RuntimeFlags) const;
    void set_flag(BodyConfigs::RuntimeFlags, bool);
    bool is_alive() const;
    void set_alive(bool);
    // Velocity of the body at a world-space point
    // (linear + angular contribution).
    Vec2 velocity_at_world_point(const Vec2&) const;
    // Linearly blend pose and velocity toward another state by fraction.
    // The caller supplies the interpolation range; no clamping is performed.
    BodyState interpolated_to(const BodyState&, real fraction) const;
};

// Mass and bounds derived from the body's attached shape geometry.
struct BodyProperties {
    real mass = 0, inv_mass = 0, inertia = 0, inv_inertia = 0;
    // Center of mass in body-local space.
    Vec2 local_center;
    // Radius enclosing all attached shapes.
    real bound_radius = 0;
    real density = 0;
    // Circle radius for lite circle bodies; box half-size for AABB bodies.
    real radius = 0;
    Vec2 half_extents;
};

// Generation-checked body identity. A reused slot receives a new generation.
struct BodyHandle {
    // Changes when a slot is reused.
    u32 generation = 0;

    constexpr BodyHandle() = default;
    constexpr BodyHandle(u32 index, u32 generation_value, BodyConfigs::Kind body_kind)
        : generation(index <= INDEX_MASK && generation_value != 0 &&
                             static_cast<u32>(body_kind) != 0 &&
                             static_cast<u32>(body_kind) <= KIND_MASK
                         ? generation_value
                         : 0),
          packed_(generation != 0 ? (static_cast<u32>(body_kind) << KIND_SHIFT) | index : 0) {}

    // Slot within the storage for kind(); use BodyPool::dense_index for column access.
    constexpr u32 index() const {
        return packed_ & INDEX_MASK;
    }
    constexpr BodyConfigs::Kind kind() const {
        return static_cast<BodyConfigs::Kind>(packed_ >> KIND_SHIFT);
    }
    bool operator==(const BodyHandle&) const = default;

private:
    static constexpr u32 KIND_SHIFT = 28;
    static constexpr u32 KIND_MASK = 0xFu;
    static constexpr u32 INDEX_MASK = (u32{1} << KIND_SHIFT) - 1;
    // High four bits contain BodyConfigs::Kind; None (zero) denotes an invalid handle.
    u32 packed_ = 0;
};
static_assert(sizeof(BodyHandle) == sizeof(u64));

class BodyPool;

// Cheap structural check; use BodyPool::is_valid for liveness.
inline constexpr bool is_valid(BodyHandle h) {
    return h.generation != 0 && h.kind() != BodyConfigs::Kind::None;
}

struct IslandId {
    // Generation guards against a reused island slot.
    u32 index = std::numeric_limits<u32>::max();
    u32 generation = 0;

    constexpr bool valid() const {
        return index != std::numeric_limits<u32>::max();
    }

    constexpr bool operator==(const IslandId& other) const = default;
};

struct Island {
    // Last simulation time reached by every body in this island.
    tick_t time = 0;

    struct SmallIsland {
        std::array<BodyHandle, 3> data;
    };

    std::variant<SmallIsland, std::vector<BodyHandle>> bodies;
};

} // namespace sqp
