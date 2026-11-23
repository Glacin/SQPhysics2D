// SPDX-FileCopyrightText: 2023 Erin Catto
// SPDX-License-Identifier: MIT

// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

// This files contains code from https://github.com/erincatto/box2d/blob/main/src/contact_solver.h

#pragma once

#include <span>
#include <type_traits>
#include <vector>

#include "sqphysics2d/body.h"
#include "sqphysics2d/detail/contact.h"

namespace sqp {
struct ContactSoftness {
    real bias_rate = 0, mass_scale = 0, impulse_scale = 0;
};

ContactSoftness make_contact_softness(real hertz, real contact_damping_ratio, real h);

struct ContactConstraintPoint {
    Vec2 anchor_a{}, anchor_b{};
    real base_separation = 0, current_separation = 0, normal_mass = 0, tangent_mass = 0;
    real normal_impulse = 0, tangent_impulse = 0;
    real restitution_velocity = 0, initial_normal_velocity = 0;
    real total_normal_impulse = 0, total_tangent_impulse = 0, peak_normal_impulse = 0;
};

struct ContactConstraint {
    BodyHandle a, b;
    Vec2 normal{1, 0}, initial_delta{};
    // TOI constraints use a pair-local contact point and body positions relative to this origin.
    Vec2 local_origin{};
    bool local_geometry = false;
    real rotation_a = 0, rotation_b = 0;
    Vec2 inverse_mass_a{}, inverse_mass_b{};
    real inverse_inertia_a = 0, inverse_inertia_b = 0;
    real static_friction = 0, dynamic_friction = 0, restitution = 0;
    bool restitution_applied = false;
    real max_penetration_correction_speed = real(3), inverse_h = 0;
    ContactSoftness softness;
    i32 count = 0;
    ContactConstraintPoint points[2];
};

template <class T> struct ContactBodyState {
    T *position_x = nullptr, *position_y = nullptr, *rotation = nullptr;
    T *velocity_x = nullptr, *velocity_y = nullptr, *angular_velocity = nullptr;
    std::conditional_t<std::is_const_v<T>, const BodyConfigs::RuntimeFlags,
                       BodyConfigs::RuntimeFlags>* runtime_flags = nullptr;
    real inverse_mass_x = 0, inverse_mass_y = 0, inverse_inertia = 0;
};

struct PreparedContactConstraint {
    ContactConstraint* constraint = nullptr;
    // Bindings stay solve-local and must be rebuilt after constraint storage changes.
    ContactBodyState<real> a, b;
};

// Final contact parameters produced before constraint preparation.
struct ContactSolverParameters {
    real restitution = 0;
    real dynamic_friction = 0;
    real static_friction = 0;
    real damping_ratio = 0;
    real contact_frequency = 0;
    real max_penetration_correction_speed = 0;
};

ContactConstraint prepare_contact_constraint(const BodyPool&, BodyHandle, BodyHandle,
                                             const ContactManifold&, const ContactSolverParameters&,
                                             real h, bool persisted,
                                             real restitution_threshold = real(0.1));

// Build the world microstep constraint without subtracting quantized world contact points.
ContactConstraint prepare_toi_local_contact_constraint(const BodyPool&, BodyHandle, BodyHandle,
                                                       const ContactManifold& local_manifold,
                                                       Vec2 origin, const ContactSolverParameters&,
                                                       real h, bool persisted,
                                                       real restitution_threshold = real(0.1));

// Bind a prepared constraint to live body columns; rebuild after storage changes.
bool prepare_contact_runtime(BodyPool&, ContactConstraint&, PreparedContactConstraint&);
void refresh_contact_constraint_geometry(PreparedContactConstraint&);
void warm_start_contact(PreparedContactConstraint&);

void solve_contact_constraint(PreparedContactConstraint&, real slop, bool use_bias);

void apply_contact_restitution(PreparedContactConstraint&);

ContactSolveReport report_contact_constraint(const ContactConstraint&);

struct ToiContact {
    BodyHandle a, b;
    ContactManifold manifold;
    ContactSolverParameters parameters;
};

// Standalone TOI solve path; World microsteps use prepare_toi_local_contact_constraint instead.
// Reports stay in input order, including contacts with no movable response.
bool solve_toi_contacts(BodyPool&, std::span<const ToiContact>, u32 velocity_iterations,
                        real restitution_threshold, std::vector<ContactSolveReport>&);

} // namespace sqp
