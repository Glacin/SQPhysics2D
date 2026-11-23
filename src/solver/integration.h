// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#pragma once

#include "sqphysics2d/detail/math.h"

namespace sqp {

// Predicts the state of a body after dt simulation time units, using the current state and
// acceleration. Semi-implicit Euler integration is used for velocity and position updates.

// Velocity Verlex might be better since interval is not in fixed frequency.
void integrate_body_velocity(Vec2& velocity, real& angular_velocity, const Vec2 acceleration,
                             const real angular_acceleration, const tick_t dt);
void integrate_body_velocity_seconds(Vec2& velocity, real& angular_velocity,
                                     const Vec2 acceleration, const real angular_acceleration,
                                     real dt_seconds);

void integrate_body_position(Vec2& position, real& angle, const Vec2& velocity,
                             const real angular_velocity, const tick_t dt);
void integrate_body_position_seconds(Vec2& position, real& angle, const Vec2& velocity,
                                     const real angular_velocity, real dt_seconds);

void integrate_body_state(Vec2& position, real& angle, Vec2& velocity, real& angular_velocity,
                          const Vec2 acceleration, const real angular_acceleration,
                          const tick_t dt);

void predict_body_state(const Vec2& start_position, const real start_angle,
                        const Vec2& start_velocity, const real start_angular_velocity,
                        const Vec2& acceleration, const real angular_acceleration, const tick_t dt,
                        Vec2& out_position, real& out_angle, Vec2& out_velocity,
                        real& out_angular_velocity);

void predict_body_state(const Vec2& start_position, const real start_angle,
                        const Vec2& start_velocity, const real start_angular_velocity,
                        const Vec2& acceleration, const real angular_acceleration, const tick_t dt,
                        Vec2& out_position, real& out_angle);

} // namespace sqp
