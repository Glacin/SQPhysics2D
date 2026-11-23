// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cmath>

#include "sqphysics2d/detail/math.h"
#include "sqphysics2d/material.h"

// TODO: after publish find a better priority system.
namespace sqp::detail {

struct CombinedContactMaterial {
    real contact_frequency_hz = real(64), max_penetration_correction_speed = real(3);
    real contact_damping_ratio = real(8);
    real restitution = real(0), restitution_impact_threshold = real(1);
    real restitution_impact_exponent = real(1);
    real dynamic_friction = real(0), static_friction = real(0);
    u8 restitution_flags = 0;
};

inline CombinedContactMaterial combine_contact_materials(const PhysicsMaterial& a,
                                                         const PhysicsMaterial& b) {
    CombinedContactMaterial m;
    if (a.friction_priority == b.friction_priority) {
        m.dynamic_friction = std::min(a.dynamic_friction, b.dynamic_friction);
        m.static_friction = std::min(a.static_friction, b.static_friction);
    } else {
        const auto& selected = a.friction_priority > b.friction_priority ? a : b;
        m.dynamic_friction = selected.dynamic_friction;
        m.static_friction = selected.static_friction;
    }
    m.static_friction = std::max(m.static_friction, m.dynamic_friction);

    if (a.restitution_priority == b.restitution_priority) {
        m.restitution = std::min(a.base_restitution, b.base_restitution);
        m.restitution_impact_threshold =
            std::min(a.restitution_impact_threshold, b.restitution_impact_threshold);
        m.restitution_impact_exponent =
            std::max(a.restitution_impact_exponent, b.restitution_impact_exponent);
        m.restitution_flags = static_cast<u8>(a.restitution_flags | b.restitution_flags);
    } else {
        const auto& selected = a.restitution_priority > b.restitution_priority ? a : b;
        m.restitution = selected.base_restitution;
        m.restitution_impact_threshold = selected.restitution_impact_threshold;
        m.restitution_impact_exponent = selected.restitution_impact_exponent;
        m.restitution_flags = selected.restitution_flags;
    }

    if (a.contact_response_priority == b.contact_response_priority) {
        m.contact_frequency_hz = std::min(a.contact_frequency_hz, b.contact_frequency_hz);
        m.max_penetration_correction_speed =
            std::min(a.max_penetration_correction_speed, b.max_penetration_correction_speed);
    } else {
        const auto& selected = a.contact_response_priority > b.contact_response_priority ? a : b;
        m.contact_frequency_hz = selected.contact_frequency_hz;
        m.max_penetration_correction_speed = selected.max_penetration_correction_speed;
    }
    if (a.contact_damping_priority == b.contact_damping_priority)
        m.contact_damping_ratio =
            a.contact_damping_ratio / real(2) + b.contact_damping_ratio / real(2);
    else {
        const auto& selected = a.contact_damping_priority > b.contact_damping_priority ? a : b;
        m.contact_damping_ratio = selected.contact_damping_ratio;
    }
    return m;
}

// This model tries to make these rigid bodies more like made by actual metal.

inline real contact_restitution(const CombinedContactMaterial& material, real closing,
                                real inverse_mass, real restitution_impact_scale, bool enabled) {
    real restitution = material.restitution;
    if (enabled &&
        !math::has_any_flags(material.restitution_flags, RESTITUTION_DISABLE_IMPACT_SCALING)) {
        const real impact =
            inverse_mass > real(0) ? closing / inverse_mass * restitution_impact_scale : real(0);

        // No impact is available to scale. Avoid division by zero and square roots of
        // negative values for stationary or separating contacts.
        if (!(impact > real(0)))
            return real(0);

        if (material.restitution_impact_exponent == real(1.0)) {
            restitution *= material.restitution_impact_threshold / impact;
        } else if (material.restitution_impact_exponent == real(0.5)) {
            restitution *= math::sqrt(material.restitution_impact_threshold / impact);
        } else if (material.restitution_impact_exponent == real(0.25)) {
            restitution *= math::sqrt(math::sqrt(material.restitution_impact_threshold / impact));
        } else if (closing > real(0) && impact > material.restitution_impact_threshold)
            restitution *= std::pow(material.restitution_impact_threshold / impact,
                                    material.restitution_impact_exponent);
    }
    return std::clamp(restitution, real(0), real(1));
}

} // namespace sqp::detail
