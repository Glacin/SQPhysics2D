// SPDX-FileCopyrightText: 2026 Erin Catto
// SPDX-License-Identifier: MIT

// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

// This files contains code from https://github.com/erincatto/box2d/blob/main/src/contact_solver.c

// Fixed-anchor contact preparation and normal push/relax are adapted from the
// vendored Box2D contact solver. SQP retains separate static and dynamic friction.

#include "sqphysics2d/detail/contact_solver.h"
#include "sqphysics2d/detail/data_storage.h"

#include <algorithm>
#include <cmath>

#include "sqphysics2d/detail/math.h"

namespace sqp {

namespace {

template <class Pool, class T>
bool load_body(Pool& bodies, BodyHandle handle, ContactBodyState<T>& out) {
    out = {};
    const auto index = bodies.dense_index(handle);
    if (!index)
        return false;
    const auto field = [&](BodyScalarColumn column) -> T* {
        auto values = bodies.scalar_column(handle.kind(), column);
        return values.empty() ? nullptr : &values[*index];
    };
    out.position_x = field(BodyScalarColumn::PositionX);
    out.position_y = field(BodyScalarColumn::PositionY);
    out.velocity_x = field(BodyScalarColumn::VelocityX);
    out.velocity_y = field(BodyScalarColumn::VelocityY);
    out.rotation = field(BodyScalarColumn::Rotation);
    out.angular_velocity = field(BodyScalarColumn::AngularVelocity);
    auto& configs = bodies.configs_column(handle.kind())[*index];
    out.runtime_flags = &configs.runtime_flags;
    const auto overrides = configs.override_flags;
    const real inverse_mass = *field(BodyScalarColumn::InverseMass);
    const auto locked = [overrides](BodyConfigs::OverrideFlags flag) {
        return (overrides & flag) == flag;
    };
    out.inverse_mass_x = locked(BodyConfigs::OverrideFlags::LockLinearX) ? real(0) : inverse_mass;
    out.inverse_mass_y = locked(BodyConfigs::OverrideFlags::LockLinearY) ? real(0) : inverse_mass;
    if (T* inertia = field(BodyScalarColumn::InverseInertia))
        out.inverse_inertia = locked(BodyConfigs::OverrideFlags::LockRot) ? real(0) : *inertia;
    return true;
}

inline real rotation(const ContactBodyState<const real>& body) {
    return body.rotation ? *body.rotation : real(0);
}
template <class T> inline real angular_velocity(const ContactBodyState<T>& body) {
    return body.angular_velocity ? *body.angular_velocity : real(0);
}
void rotate(real x, real y, real angle, real& out_x, real& out_y) {
    const real c = std::cos(angle), s = std::sin(angle);
    out_x = c * x - s * y;
    out_y = s * x + c * y;
}

Vec2 rotated(Vec2 value, real angle) {
    Vec2 result;
    rotate(value.x, value.y, angle, result.x, result.y);
    return result;
}

Vec2 contact_body_delta(const PreparedContactConstraint& prepared) {
    const auto& c = *prepared.constraint;
    if (c.local_geometry)
        return {(*prepared.b.position_x - c.local_origin.x) -
                    (*prepared.a.position_x - c.local_origin.x),
                (*prepared.b.position_y - c.local_origin.y) -
                    (*prepared.a.position_y - c.local_origin.y)};
    return {*prepared.b.position_x - *prepared.a.position_x,
            *prepared.b.position_y - *prepared.a.position_y};
}

real inverse_effective_mass(real axis_x, real axis_y, real mass_x, real mass_y, real inertia_a,
                            real inertia_b, real rax, real ray, real rbx, real rby) {
    const real ca = rax * axis_y - ray * axis_x;
    const real cb = rbx * axis_y - rby * axis_x;
    const real k = axis_x * mass_x * axis_x + axis_y * mass_y * axis_y + inertia_a * ca * ca +
                   inertia_b * cb * cb;
    return k > config::EPSILON * config::EPSILON ? real(1) / k : real(0);
}

template <class A, class B>
real normal_velocity(const A& a, const B& b, real rax, real ray, real rbx, real rby, real nx,
                     real ny) {
    const real av = angular_velocity(a), bv = angular_velocity(b);
    const real relative_x = *b.velocity_x - bv * rby - *a.velocity_x + av * ray;
    const real relative_y = *b.velocity_y + bv * rbx - *a.velocity_y - av * rax;
    return relative_x * nx + relative_y * ny;
}

struct ContactVelocity {
    Vec2 linear;
    real angular;
    Vec2 inverse_mass;
    real inverse_inertia;
};

ContactVelocity load_velocity(const ContactBodyState<real>& body) {
    return {{*body.velocity_x, *body.velocity_y},
            body.angular_velocity ? *body.angular_velocity : real(0),
            {body.inverse_mass_x, body.inverse_mass_y},
            body.inverse_inertia};
}

void store_velocity(const ContactBodyState<real>& body, const ContactVelocity& velocity) {
    *body.velocity_x = velocity.linear.x;
    *body.velocity_y = velocity.linear.y;
    if (body.angular_velocity)
        *body.angular_velocity = velocity.angular;
}

real normal_velocity(const ContactVelocity& a, const ContactVelocity& b, real rax, real ray,
                     real rbx, real rby, real nx, real ny) {
    return (b.linear.x - b.angular * rby - a.linear.x + a.angular * ray) * nx +
           (b.linear.y + b.angular * rbx - a.linear.y - a.angular * rax) * ny;
}

void apply_impulse(ContactVelocity& a, ContactVelocity& b, real rax, real ray, real rbx, real rby,
                   real impulse_x, real impulse_y) {
    a.linear.x -= a.inverse_mass.x * impulse_x;
    a.linear.y -= a.inverse_mass.y * impulse_y;
    b.linear.x += b.inverse_mass.x * impulse_x;
    b.linear.y += b.inverse_mass.y * impulse_y;
    a.angular -= a.inverse_inertia * (rax * impulse_y - ray * impulse_x);
    b.angular += b.inverse_inertia * (rbx * impulse_y - rby * impulse_x);
}

} // namespace

ContactSoftness make_contact_softness(real hertz, real contact_damping_ratio, real h) {
    if (!(hertz > real(0)) || !(h > real(0)))
        return {};
    const real omega = real(2) * math::PI * hertz;
    const real a1 = real(2) * contact_damping_ratio + h * omega;
    const real a2 = h * omega * a1;
    if (!std::isfinite(omega) || !std::isfinite(a1) || !std::isfinite(a2))
        return {real(1) / h, real(1), real(0)};
    const real a3 = real(1) / (real(1) + a2);
    return {omega / a1, a2 * a3, a3};
}

static ContactConstraint prepare_contact_constraint_impl(
    const BodyPool& bodies, BodyHandle a_handle, BodyHandle b_handle,
    const ContactManifold& manifold, const ContactSolverParameters& parameters, real h,
    bool persisted, real restitution_threshold, const Vec2* local_origin) {
    ContactConstraint c;
    c.a = a_handle;
    c.b = b_handle;
    c.local_geometry = local_origin != nullptr;
    if (local_origin)
        c.local_origin = *local_origin;
    ContactBodyState<const real> a, b;
    constexpr auto alive = BodyConfigs::RuntimeFlags::IsAlive;
    if (!(h > real(0)) || manifold.count <= 0 || !load_body(bodies, a_handle, a) ||
        !load_body(bodies, b_handle, b) ||
        ((*a.runtime_flags & alive) == BodyConfigs::RuntimeFlags(0)) ||
        ((*b.runtime_flags & alive) == BodyConfigs::RuntimeFlags(0)))
        return c;

    const real normal_len2 =
        manifold.normal.x * manifold.normal.x + manifold.normal.y * manifold.normal.y;
    if (normal_len2 > config::EPSILON * config::EPSILON) {
        const real inv_len = real(1) / std::sqrt(normal_len2);
        c.normal.x = manifold.normal.x * inv_len;
        c.normal.y = manifold.normal.y * inv_len;
    }
    const Vec2 position_a = local_origin
                                ? Vec2{*a.position_x - local_origin->x,
                                       *a.position_y - local_origin->y}
                                : Vec2{*a.position_x, *a.position_y};
    const Vec2 position_b = local_origin
                                ? Vec2{*b.position_x - local_origin->x,
                                       *b.position_y - local_origin->y}
                                : Vec2{*b.position_x, *b.position_y};
    c.initial_delta = position_b - position_a;
    c.rotation_a = rotation(a);
    c.rotation_b = rotation(b);
    c.inverse_mass_a.x = a.inverse_mass_x;
    c.inverse_mass_a.y = a.inverse_mass_y;
    c.inverse_mass_b.x = b.inverse_mass_x;
    c.inverse_mass_b.y = b.inverse_mass_y;
    c.inverse_inertia_a = a.inverse_inertia;
    c.inverse_inertia_b = b.inverse_inertia;

    c.static_friction = parameters.static_friction;
    c.dynamic_friction = parameters.dynamic_friction;
    c.restitution = parameters.restitution;
    c.max_penetration_correction_speed = parameters.max_penetration_correction_speed;
    c.inverse_h = real(1) / h;
    const real hz = parameters.contact_frequency;
    c.softness = make_contact_softness((std::max(a.inverse_mass_x, a.inverse_mass_y) == real(0) ||
                                        std::max(b.inverse_mass_x, b.inverse_mass_y) == real(0))
                                           ? real(2) * hz
                                           : hz,
                                       parameters.damping_ratio, h);
    c.count = std::min(manifold.count, i32(2));
    const real tangent_x = -c.normal.y, tangent_y = c.normal.x;
    for (i32 i = 0; i < c.count; ++i) {
        ContactConstraintPoint& p = c.points[i];
        p.anchor_a = manifold.points[i].point - position_a;
        p.anchor_b = manifold.points[i].point - position_b;
        p.base_separation = -manifold.points[i].depth;
        p.current_separation = p.base_separation;
        p.normal_mass = inverse_effective_mass(
            c.normal.x, c.normal.y, c.inverse_mass_a.x + c.inverse_mass_b.x,
            c.inverse_mass_a.y + c.inverse_mass_b.y, c.inverse_inertia_a, c.inverse_inertia_b,
            p.anchor_a.x, p.anchor_a.y, p.anchor_b.x, p.anchor_b.y);
        p.tangent_mass = inverse_effective_mass(
            tangent_x, tangent_y, c.inverse_mass_a.x + c.inverse_mass_b.x,
            c.inverse_mass_a.y + c.inverse_mass_b.y, c.inverse_inertia_a, c.inverse_inertia_b,
            p.anchor_a.x, p.anchor_a.y, p.anchor_b.x, p.anchor_b.y);
        p.initial_normal_velocity = normal_velocity(a, b, p.anchor_a.x, p.anchor_a.y, p.anchor_b.x,
                                                    p.anchor_b.y, c.normal.x, c.normal.y);
        // A speculative point with a positive gap may need a stopping impulse, but it has not
        // reached an impact pose and must not receive restitution yet.
        if (!persisted && p.base_separation <= real(0) &&
            p.initial_normal_velocity < -restitution_threshold)
            p.restitution_velocity = -c.restitution * p.initial_normal_velocity;
    }
    return c;
}

ContactConstraint prepare_contact_constraint(const BodyPool& bodies, BodyHandle a_handle,
                                             BodyHandle b_handle, const ContactManifold& manifold,
                                             const ContactSolverParameters& parameters, real h,
                                             bool persisted, real restitution_threshold) {
    return prepare_contact_constraint_impl(bodies, a_handle, b_handle, manifold, parameters, h,
                                           persisted, restitution_threshold, nullptr);
}

ContactConstraint prepare_toi_local_contact_constraint(
    const BodyPool& bodies, BodyHandle a_handle, BodyHandle b_handle,
    const ContactManifold& local_manifold, Vec2 origin,
    const ContactSolverParameters& parameters, real h, bool persisted,
    real restitution_threshold) {
    return prepare_contact_constraint_impl(bodies, a_handle, b_handle, local_manifold, parameters,
                                           h, persisted, restitution_threshold, &origin);
}

bool prepare_contact_runtime(BodyPool& bodies, ContactConstraint& c,
                             PreparedContactConstraint& out) {
    ContactBodyState<real> a, b;
    if (!load_body(bodies, c.a, a) || !load_body(bodies, c.b, b)) {
        out = {};
        return false;
    }
    out = {&c, a, b};
    return true;
}

void refresh_contact_constraint_geometry(PreparedContactConstraint& prepared) {
    if (!prepared.constraint)
        return;
    auto& c = *prepared.constraint;
    const auto& a = prepared.a;
    const auto& b = prepared.b;
    const real angle_a = (a.rotation ? *a.rotation : real(0)) - c.rotation_a;
    const real angle_b = (b.rotation ? *b.rotation : real(0)) - c.rotation_b;
    const Vec2 delta = contact_body_delta(prepared);
    const Vec2 translation = delta - c.initial_delta;
    for (i32 i = 0; i < c.count; ++i) {
        auto& p = c.points[i];
        const Vec2 moved_a = rotated(p.anchor_a, angle_a);
        const Vec2 moved_b = rotated(p.anchor_b, angle_b);
        const Vec2 displacement = translation + moved_b - p.anchor_b - moved_a + p.anchor_a;
        p.current_separation = p.base_separation +
                               displacement.x * c.normal.x + displacement.y * c.normal.y;
    }
}

void warm_start_contact(PreparedContactConstraint& prepared) {
    if (!prepared.constraint)
        return;
    auto a = load_velocity(prepared.a), b = load_velocity(prepared.b);
    auto& c = *prepared.constraint;
    const real tx = -c.normal.y, ty = c.normal.x;
    for (i32 i = 0; i < c.count; ++i) {
        ContactConstraintPoint& p = c.points[i];
        p.peak_normal_impulse = std::max(p.peak_normal_impulse, p.normal_impulse);
        p.total_normal_impulse += p.normal_impulse;
        p.total_tangent_impulse += std::abs(p.tangent_impulse);
        apply_impulse(a, b, p.anchor_a.x, p.anchor_a.y, p.anchor_b.x, p.anchor_b.y,
                      c.normal.x * p.normal_impulse + tx * p.tangent_impulse,
                      c.normal.y * p.normal_impulse + ty * p.tangent_impulse);
    }
    store_velocity(prepared.a, a);
    store_velocity(prepared.b, b);
}

void solve_contact_constraint(PreparedContactConstraint& prepared, real slop, bool use_bias) {
    if (!(prepared.constraint && prepared.constraint->inverse_h > real(0)) ||
        prepared.constraint->count <= 0)
        return;
    auto a = load_velocity(prepared.a), b = load_velocity(prepared.b);
    auto& c = *prepared.constraint;
    const real tx = -c.normal.y, ty = c.normal.x;
    real normal_impulses[2]{};
    real normal_cross_a[2]{}, normal_cross_b[2]{};
    real normal_impulse = 0, angular_impulse_a = 0, angular_impulse_b = 0;
    const real normal_linear_mass =
        (c.inverse_mass_a.x + c.inverse_mass_b.x) * c.normal.x * c.normal.x +
        (c.inverse_mass_a.y + c.inverse_mass_b.y) * c.normal.y * c.normal.y;
    const real base_normal_velocity =
        (b.linear.x - a.linear.x) * c.normal.x +
        (b.linear.y - a.linear.y) * c.normal.y;
    for (i32 i = 0; i < c.count; ++i) {
        ContactConstraintPoint& p = c.points[i];
        const real separation = p.current_separation;
        real bias = 0, mass_scale = 1, impulse_scale = 0;
        if (separation > real(0))
            bias = separation * c.inverse_h;
        else if (use_bias) {
            mass_scale = c.softness.mass_scale;
            impulse_scale = c.softness.impulse_scale;
            bias =
                std::max(mass_scale * c.softness.bias_rate * std::min(real(0), separation + slop),
                         -c.max_penetration_correction_speed);
        }
        normal_cross_a[i] = p.anchor_a.x * c.normal.y - p.anchor_a.y * c.normal.x;
        normal_cross_b[i] = p.anchor_b.x * c.normal.y - p.anchor_b.y * c.normal.x;
        real vn = base_normal_velocity + b.angular * normal_cross_b[i] -
                  a.angular * normal_cross_a[i];
        for (i32 j = 0; j < i; ++j)
            vn += (normal_linear_mass +
                   c.inverse_inertia_a * normal_cross_a[i] * normal_cross_a[j] +
                   c.inverse_inertia_b * normal_cross_b[i] * normal_cross_b[j]) *
                  normal_impulses[j];
        const real previous = p.normal_impulse;
        p.normal_impulse = std::max(real(0), previous - p.normal_mass * (mass_scale * vn + bias) -
                                                 impulse_scale * previous);
        p.peak_normal_impulse = std::max(p.peak_normal_impulse, p.normal_impulse);
        const real impulse = p.normal_impulse - previous;
        p.total_normal_impulse += impulse;
        normal_impulses[i] = impulse;
        normal_impulse += impulse;
        angular_impulse_a += impulse * normal_cross_a[i];
        angular_impulse_b += impulse * normal_cross_b[i];
    }
    apply_impulse(a, b, 0, 0, 0, 0,
                  c.normal.x * normal_impulse, c.normal.y * normal_impulse);
    a.angular -= c.inverse_inertia_a * angular_impulse_a;
    b.angular += c.inverse_inertia_b * angular_impulse_b;
    if (use_bias) {
        store_velocity(prepared.a, a);
        store_velocity(prepared.b, b);
        return;
    }
    real total_normal_impulse = 0;
    real total_tangent_impulse = 0;
    Vec2 average_a{}, average_b{};
    for (i32 i = 0; i < c.count; ++i) {
        const auto& point = c.points[i];
        total_normal_impulse += point.normal_impulse;
        total_tangent_impulse += point.tangent_impulse;
        average_a += point.anchor_a;
        average_b += point.anchor_b;
    }
    average_a /= static_cast<real>(c.count);
    average_b /= static_cast<real>(c.count);
    const real pair_tangent_mass = inverse_effective_mass(
        tx, ty, c.inverse_mass_a.x + c.inverse_mass_b.x,
        c.inverse_mass_a.y + c.inverse_mass_b.y, c.inverse_inertia_a, c.inverse_inertia_b,
        average_a.x, average_a.y, average_b.x, average_b.y);
    const real pair_tangent_velocity = normal_velocity(
        a, b, average_a.x, average_a.y, average_b.x, average_b.y, tx, ty);
    const bool sliding = std::abs(total_tangent_impulse -
                                  pair_tangent_mass * pair_tangent_velocity) >
                         c.static_friction * total_normal_impulse;
    for (i32 i = 0; i < c.count; ++i) {
        ContactConstraintPoint& p = c.points[i];
        const real vt = normal_velocity(a, b, p.anchor_a.x, p.anchor_a.y, p.anchor_b.x,
                                        p.anchor_b.y, tx, ty);
        const real previous = p.tangent_impulse;
        const real velocity_candidate = previous - p.tangent_mass * vt;
        const real dynamic_limit = c.dynamic_friction * p.normal_impulse;
        const real static_limit = c.static_friction * p.normal_impulse;
        if (!sliding) {
            p.tangent_impulse = std::clamp(velocity_candidate, -static_limit, static_limit);
        }
        else {
            p.tangent_impulse = std::clamp(velocity_candidate, -dynamic_limit, dynamic_limit);
        }
        const real impulse = p.tangent_impulse - previous;
        p.total_tangent_impulse += std::abs(impulse);
        apply_impulse(a, b, p.anchor_a.x, p.anchor_a.y, p.anchor_b.x, p.anchor_b.y,
                      tx * impulse, ty * impulse);
    }
    store_velocity(prepared.a, a);
    store_velocity(prepared.b, b);
}

void apply_contact_restitution(PreparedContactConstraint& prepared) {
    if (!prepared.constraint)
        return;
    auto& c = *prepared.constraint;
    if (!(c.restitution > real(0)))
        return;
    auto a = load_velocity(prepared.a), b = load_velocity(prepared.b);
    for (i32 i = 0; i < c.count; ++i) {
        ContactConstraintPoint& p = c.points[i];
        if (!(p.restitution_velocity > real(0)) || !(p.peak_normal_impulse > real(0)))
            continue;
        const real vn = normal_velocity(a, b, p.anchor_a.x, p.anchor_a.y, p.anchor_b.x,
                                        p.anchor_b.y, c.normal.x, c.normal.y);
        const real previous = p.normal_impulse;
        p.normal_impulse =
            std::max(real(0), previous - p.normal_mass * (vn - p.restitution_velocity));
        const real impulse = p.normal_impulse - previous;
        if (impulse == real(0))
            continue;
        p.peak_normal_impulse = std::max(p.peak_normal_impulse, p.normal_impulse);
        p.total_normal_impulse += impulse;
        apply_impulse(a, b, p.anchor_a.x, p.anchor_a.y, p.anchor_b.x, p.anchor_b.y,
                      c.normal.x * impulse, c.normal.y * impulse);
        c.restitution_applied = true;
    }
    store_velocity(prepared.a, a);
    store_velocity(prepared.b, b);
}

ContactSolveReport report_contact_constraint(const ContactConstraint& c) {
    ContactSolveReport report;
    report.restitution_applied = c.restitution_applied;
    for (i32 i = 0; i < c.count; ++i) {
        const ContactConstraintPoint& p = c.points[i];
        report.points[i] = {p.normal_impulse, p.tangent_impulse};
        report.normal_impulse += std::max(real(0), p.total_normal_impulse);
        report.tangent_impulse += p.total_tangent_impulse;
    }
    report.total_impulse = report.normal_impulse + report.tangent_impulse;
    return report;
}

bool solve_toi_contacts(BodyPool& bodies, std::span<const ToiContact> contacts,
                        u32 velocity_iterations, real restitution_threshold,
                        std::vector<ContactSolveReport>& reports) {
    reports.clear();
    std::vector<ContactConstraint> constraints;
    constraints.reserve(contacts.size());
    for (const auto& contact : contacts) {
        if (contact.manifold.count == 0)
            continue;
        // TOI is a velocity-only solve. The preparation interval is irrelevant to
        // softness because the bias pass and ordinary position integration are omitted.
        auto constraint =
            prepare_contact_constraint(bodies, contact.a, contact.b, contact.manifold,
                                       contact.parameters, real(1), false, restitution_threshold);
        if (constraint.count != std::clamp(contact.manifold.count, 0, 2))
            return false;
        constraints.push_back(std::move(constraint));
    }
    std::vector<PreparedContactConstraint> prepared(constraints.size());
    for (std::size_t i = 0; i < constraints.size(); ++i)
        if (!prepare_contact_runtime(bodies, constraints[i], prepared[i]))
            return false;
    for (u32 pass = 0; pass < std::max(1u, velocity_iterations); ++pass)
        for (auto& constraint : prepared)
            solve_contact_constraint(constraint, real(0), false);
    reports.reserve(prepared.size());
    for (auto& constraint : prepared) {
        apply_contact_restitution(constraint);
        reports.push_back(report_contact_constraint(*constraint.constraint));
    }
    return true;
}

} // namespace sqp
