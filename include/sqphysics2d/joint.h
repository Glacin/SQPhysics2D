// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#pragma once

#include <limits>

#include "sqphysics2d/body.h"

namespace sqp {

// Joint definitions and stored joint data share these five constraint types.
// The current World pipeline does not submit joints to the contact solver.
enum class JointType : u8 { Pin, Slider, Weld, Distance, Spring };

// A slot remains stable while the joint is alive. A reused slot gets a new generation.
struct JointHandle {
    u32 index = std::numeric_limits<u32>::max();
    u32 generation = 0;

    constexpr bool valid() const {
        return index <= static_cast<u32>(std::numeric_limits<i32>::max() >> 1) && generation != 0;
    }
    constexpr bool operator==(const JointHandle&) const = default;
    constexpr auto operator<=>(const JointHandle&) const = default;
};

// Common joint endpoints. An invalid A body denotes the fixed world origin;
// B must be a live body. Anchors are local to their respective endpoints.
struct JointDef {
    BodyHandle body_a{};
    BodyHandle body_b{};
    Vec2 local_anchor_a{};
    Vec2 local_anchor_b{};
    // Whether the connected bodies remain eligible for collision candidates.
    bool collide_connected = false;
};

// Defines a shared anchor while permitting relative rotation.
struct PinJointDef : JointDef {};

// Defines translation along A's local axis, optionally within the given limits.
struct SliderJointDef : JointDef {
    Vec2 local_axis_a{1, 0};
    bool enable_limit = false;
    real lower_translation = 0;
    real upper_translation = 0;
};

// Defines fixed anchors and relative angle.
struct WeldJointDef : JointDef {};

// Defines a fixed anchor-to-anchor distance in metres.
struct DistanceJointDef : JointDef {
    real length = 0;
};

// Defines a spring with rest_length, frequency in Hz and dimensionless damping ratio.
struct SpringJointDef : JointDef {
    real rest_length = 0;
    real hertz = 5;
    real damping_ratio = real(0.7);
};

// Internal adjacency list slot. Each live joint has up to two endpoint edges.
inline constexpr u32 INVALID_JOINT_EDGE = std::numeric_limits<u32>::max();
inline constexpr u32 max_joint_slot = static_cast<u32>(std::numeric_limits<i32>::max() >> 1);

struct JointEdge {
    u32 prev = INVALID_JOINT_EDGE;
    u32 next = INVALID_JOINT_EDGE;
};

// Endpoint zero belongs to A, endpoint one to B. The world anchor has no A edge.
constexpr u32 joint_edge_key(JointHandle handle, u32 side) {
    return handle.valid() && side < 2 ? (handle.index << 1) | side : INVALID_JOINT_EDGE;
}

struct JointRegistration {
    // Pool-owned endpoint links and index into the type-specific storage.
    BodyHandle body_a{};
    BodyHandle body_b{};
    JointEdge edges[2]{};
    u32 type_index = 0;
    JointType type = JointType::Pin;
    bool collide_connected = false;
};

// Type-specific runtime values captured when a joint is created.
struct PinJoint {
    Vec2 local_anchor_a{};
    Vec2 local_anchor_b{};
};

struct SliderJoint {
    Vec2 local_anchor_a{};
    Vec2 local_anchor_b{};
    Vec2 local_axis_a{1, 0};
    real reference_angle = 0;
    real lower_translation = 0;
    real upper_translation = 0;
    bool enable_limit = false;
};

struct WeldJoint {
    Vec2 local_anchor_a{};
    Vec2 local_anchor_b{};
    real reference_angle = 0;
};

struct DistanceJoint {
    Vec2 local_anchor_a{};
    Vec2 local_anchor_b{};
    real length = 0;
};

struct SpringJoint {
    Vec2 local_anchor_a{};
    Vec2 local_anchor_b{};
    real rest_length = 0;
    real hertz = 5;
    real damping_ratio = real(0.7);
};

// Validate live endpoints, finite anchors and each type's extra parameters.
// These functions do not modify the definition or create a joint.
bool validate_joint_def(const PinJointDef&, const BodyPool&);
bool validate_joint_def(const SliderJointDef&, const BodyPool&);
bool validate_joint_def(const WeldJointDef&, const BodyPool&);
bool validate_joint_def(const DistanceJointDef&, const BodyPool&);
bool validate_joint_def(const SpringJointDef&, const BodyPool&);

} // namespace sqp
