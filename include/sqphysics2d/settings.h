// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#pragma once

#include "sqphysics2d/detail/math.h"

namespace sqp {

// Bounds the scheduler's ordinary intervals and the number of TOI interruptions.
struct IntervalSettings {
    // Maximum predicted translation and rotation for bodies not using CCD.
    real max_linear_displacement_m = real(0.1);
    real max_angular_displacement_rad = real(0.01);
    // Minimum period for every body kind, rounded up to a power of two TU.
    tick_t min_period_tu = config::TU_PER_TICK / 4;
    // Maximum TOI truncations charged to an island during one nominal tick.
    u32 max_truncations_per_tick = 4;
    // Replace invalid values with defaults and clamp the period and budget.
    void sanitize() noexcept;
};

// Settings for the ordinary and world microstep contact solver.
struct ContactSettings {
    // Ordinary contacts use substeps; TOI microsteps use one integration step.
    u32 substeps = 4;
    // Number of velocity passes in both ordinary and world microstep solves.
    u32 velocity_iterations = 1;

    // Used by the separate TOI contact solver, not by World microsteps.
    u32 toi_velocity_iterations = 4;

    // Contacts below this closing speed do not receive a bounce impulse.
    real restitution_threshold_m_s = real(0.1);
    // Enable impact-dependent scaling of material restitution.
    bool enable_dynamic_restitution = true;
};

// Selects which shape families may use CCD and when automatic casts begin.
struct CcdSettings {
    bool enable_circle = true, enable_polygon = true;
    // Auto compares translation and angular travel at the shape's sweep radius
    // against this shared distance threshold.
    real linear_start_threshold_m = real(0.4);
};

// Extra tree-node space around each dynamic shape's exact swept AABB.
struct BroadPhaseSettings {
    // Each side grows by this fraction of the swept AABB's half-width or half-height.
    // Zero disables the extra space. The default is selected by the benchmark matrix.
    static constexpr real default_proxy_extension_ratio = real(0.2);
    real proxy_extension_ratio = default_proxy_extension_ratio;
};

// Complete settings for a World. Values are sanitized when set by command.
struct PhysicsSettings {
    IntervalSettings interval;
    ContactSettings contact;
    CcdSettings ccd;
    BroadPhaseSettings broad_phase;

    // World-space acceleration in metres per second squared.
    Vec2 gravity{0, real(-9.81)};
    // Stored and sanitized for compatibility; the current world integration does not read them.
    real linear_velocity_tick_multiplier = 1;
    real angular_velocity_tick_multiplier = 1;

    // Sleep after an island stays below both speed thresholds for the quiet duration.
    bool allow_sleep = true;
    real sleep_linear_velocity_threshold_m_s = real(0.01);
    real sleep_angular_velocity_threshold_rad_s = real(0.05);

    tick_t sleep_time_threshold_tu = config::TU_PER_TICK * 16;
    // Clamp counts and replace invalid physical values with defaults.
    void sanitize() noexcept;
};

} // namespace sqp
