// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "sqphysics2d/body.h"
#include "sqphysics2d/detail/math.h"

namespace sqp {

Vec2 gravity_with_body_overrides(Vec2 gravity, BodyConfigs::OverrideFlags flags) {
    gravity.x *= !math::has_any_flags(flags, BodyConfigs::OverrideFlags::NoGravityX);
    gravity.y *= !math::has_any_flags(flags, BodyConfigs::OverrideFlags::NoGravityY);
    return gravity;
}

RigidBodyDef::RigidBodyDef(const Vec2& p, real r, const Vec2& v, real w)
    : position(p), rotation(r), linear_velocity(v), angular_velocity(w) {}

LiteCircleBodyDef::LiteCircleBodyDef(const Vec2& p, const Vec2& v, real drag)
    : position(p), linear_velocity(v), drag_multiplier(drag) {}

AabbBoxBodyDef::AabbBoxBodyDef(const Vec2& p, const Vec2& v, real drag)
    : position(p), linear_velocity(v), drag_multiplier(drag) {}

bool BodyState::has_flag(BodyConfigs::RuntimeFlags flag) const {
    return (flags & flag) == flag;
}

void BodyState::set_flag(BodyConfigs::RuntimeFlags flag, bool enabled) {
    flags = enabled ? flags | flag : flags & ~flag;
}

bool BodyState::is_alive() const {
    return has_flag(BodyConfigs::RuntimeFlags::IsAlive);
}

void BodyState::set_alive(bool value) {
    set_flag(BodyConfigs::RuntimeFlags::IsAlive, value);
}

Vec2 BodyState::velocity_at_world_point(const Vec2& point) const {
    return linear_velocity + math::cross(angular_velocity, point - position);
}

BodyState BodyState::interpolated_to(const BodyState& end, real fraction) const {
    BodyState result = *this;
    result.position += (end.position - position) * fraction;
    result.rotation += (end.rotation - rotation) * fraction;
    result.linear_velocity += (end.linear_velocity - linear_velocity) * fraction;
    result.angular_velocity += (end.angular_velocity - angular_velocity) * fraction;
    return result;
}

} // namespace sqp
