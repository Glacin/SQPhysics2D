// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#pragma once

#include <limits>

#include "sqphysics2d/detail/types.h"

namespace sqp {

// Bits stored in PhysicsMaterial::restitution_flags.
enum PhysicsMaterialRestitutionFlags : u8 {
    // Skip impact-dependent restitution scaling after materials are combined.
    RESTITUTION_DISABLE_IMPACT_SCALING = 1u << 0
};

// Surface response used when two shapes contact.
// Each property group selects the higher-priority material. On a tie, friction,
// base restitution, impact threshold and response use the lesser value;
// impact exponent uses the greater value, damping is averaged, and flags are ORed.
struct PhysicsMaterial {

    // Inclusive bounds of each priority byte.
    static inline u8 MAX_PRIORITY = 255;
    static inline u8 MIN_PRIORITY = 0;

    // Soft-contact frequency in Hz and maximum correction speed in m/s.
    real contact_frequency_hz = real(64);
    real max_penetration_correction_speed = real(3);
    // Dimensionless damping ratio.
    real contact_damping_ratio = real(8);
    // Base bounce coefficient and impact-dependent scaling parameters.
    real base_restitution = real(1);
    real restitution_impact_threshold = real(1);
    real restitution_impact_exponent = real(1);
    // Dimensionless friction coefficients; static should be at least dynamic.
    real dynamic_friction = real(0.5), static_friction = real(0.75);
    // Independent priorities for the four material property groups.
    u8 restitution_priority = 0, friction_priority = 0, contact_damping_priority = 0,
       contact_response_priority = 0;
    // PhysicsMaterialRestitutionFlags bits.
    u8 restitution_flags = 0;
};

// Sentinel for "no material".
inline constexpr u32 INVALID_MATERIAL_INDEX = std::numeric_limits<u32>::max();

class MaterialPool;

} // namespace sqp
