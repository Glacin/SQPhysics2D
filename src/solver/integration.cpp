// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "sqphysics2d/detail/math.h"

#include "integration.h"

namespace sqp {

void integrate_body_velocity(Vec2& velocity, real& angular_velocity, const Vec2 acceleration,
                             const real angular_acceleration, const tick_t dt) {
    integrate_body_velocity_seconds(velocity, angular_velocity, acceleration, angular_acceleration,
                                    math::tu_to_seconds<real>(dt));
}

void integrate_body_velocity_seconds(Vec2& velocity, real& angular_velocity,
                                     const Vec2 acceleration, const real angular_acceleration,
                                     const real dt_seconds) {
    velocity = velocity + acceleration * dt_seconds;
    angular_velocity = angular_velocity + angular_acceleration * dt_seconds;
}

void integrate_body_position(Vec2& position, real& angle, const Vec2& velocity,
                             const real angular_velocity, const tick_t dt) {
    integrate_body_position_seconds(position, angle, velocity, angular_velocity,
                                    math::tu_to_seconds<real>(dt));
}

void integrate_body_position_seconds(Vec2& position, real& angle, const Vec2& velocity,
                                     const real angular_velocity, const real dt_seconds) {
    position = position + velocity * dt_seconds;
    angle = angle + angular_velocity * dt_seconds;
}

void integrate_body_state(Vec2& position, real& angle, Vec2& velocity, real& angular_velocity,
                          const Vec2 acceleration, const real angular_acceleration,
                          const tick_t dt) {
    integrate_body_velocity(velocity, angular_velocity, acceleration, angular_acceleration, dt);
    integrate_body_position(position, angle, velocity, angular_velocity, dt);
}

void predict_body_state(const Vec2& start_position, const real start_angle,
                        const Vec2& start_velocity, const real start_angular_velocity,
                        const Vec2& acceleration, const real angular_acceleration, const tick_t dt,
                        Vec2& out_position, real& out_angle, Vec2& out_velocity,
                        real& out_angular_velocity) {

    out_position = start_position;
    out_angle = start_angle;
    out_velocity = start_velocity;
    out_angular_velocity = start_angular_velocity;
    integrate_body_state(out_position, out_angle, out_velocity, out_angular_velocity, acceleration,
                         angular_acceleration, dt);
};

void predict_body_state(const Vec2& start_position, const real start_angle,
                        const Vec2& start_velocity, const real start_angular_velocity,
                        const Vec2& acceleration, const real angular_acceleration, const tick_t dt,
                        Vec2& out_position, real& out_angle) {
    out_position = start_position;
    out_angle = start_angle;
    Vec2 velocity = start_velocity;
    real angular_velocity = start_angular_velocity;
    integrate_body_state(out_position, out_angle, velocity, angular_velocity, acceleration,
                         angular_acceleration, dt);
};
} // namespace sqp
