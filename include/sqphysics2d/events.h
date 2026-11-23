// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#pragma once

#include "sqphysics2d/detail/contact.h"
#include "sqphysics2d/detail/types.h"

namespace sqp {

// Contact callbacks currently emitted by the world pipeline.
enum class PhysicsEventTypes : u8 {
    PreSolve = 0,
    PostSolve,
};

// Mutable contact and combined material values passed before a pair is solved.
// A callback may set manifold.count to zero to disable this contact for the solve.
// Contact points are in world coordinates, including for TOI microsteps.
struct PreSolveEvent {
    Pair contact_pair;
    ContactManifold manifold;

    // Restitution and friction after combining the two shape materials.
    real combined_restitution = real(0);
    real combined_dynamic_friction = real(0);
    real combined_static_friction = real(0);

    // Contact softness and the maximum penetration correction speed.
    real combined_damping_ratio = real(0);
    real combined_contact_frequency = real(0);
    real combined_max_penetration_correction_speed = real(0);
};

// Solved impulse report with the original world-space manifold from before PreSolve edits.
struct PostSolveEvent {
    Pair contact_pair;
    ContactManifold manifold;
    ContactSolveReport report;
};

// Tagged event value delivered to EventCallback.
struct PhysicsEvent {
    PhysicsEventTypes type;
    std::variant<PreSolveEvent, PostSolveEvent> event;
};

// Called synchronously while stepping. user_data is the pointer registered with World.
// PreSolve edits are validated and used by the solver; PostSolve edits are not applied.
using EventCallback = void (*)(PhysicsEvent& event, void* user_data);

} // namespace sqp
