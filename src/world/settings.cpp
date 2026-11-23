// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "sqphysics2d/settings.h"

#include <algorithm>
#include <bit>
#include <cmath>

namespace sqp {

void IntervalSettings::sanitize() noexcept {
    if (!std::isfinite(max_linear_displacement_m) || max_linear_displacement_m <= 0)
        max_linear_displacement_m = real(0.1);
    if (!std::isfinite(max_angular_displacement_rad) || max_angular_displacement_rad <= 0)
        max_angular_displacement_rad = real(0.01);
    const tick_t period = min_period_tu ? min_period_tu : config::TU_PER_TICK / 4;
    min_period_tu = std::bit_ceil(std::clamp<tick_t>(period, 1, config::TU_PER_TICK));
    max_truncations_per_tick = std::max(1u, max_truncations_per_tick);
}

void PhysicsSettings::sanitize() noexcept {
    interval.sanitize();
    const auto non_negative = [](real& value, real fallback) {
        if (!std::isfinite(value) || value < 0)
            value = fallback;
    };
    contact.substeps = std::clamp(contact.substeps, 1u, 64u);
    contact.velocity_iterations = std::clamp(contact.velocity_iterations, 1u, 256u);
    contact.toi_velocity_iterations = std::clamp(contact.toi_velocity_iterations, 1u, 256u);

    non_negative(contact.restitution_threshold_m_s, real(0.1));
    non_negative(ccd.linear_start_threshold_m, real(0.4));
    if (!std::isfinite(broad_phase.proxy_extension_ratio))
        broad_phase.proxy_extension_ratio = BroadPhaseSettings::default_proxy_extension_ratio;
    else
        broad_phase.proxy_extension_ratio =
            std::clamp(broad_phase.proxy_extension_ratio, real(0), real(1));
    if (!std::isfinite(gravity.x))
        gravity.x = 0;
    if (!std::isfinite(gravity.y))
        gravity.y = 0;

    non_negative(linear_velocity_tick_multiplier, real(1));
    non_negative(angular_velocity_tick_multiplier, real(1));
    non_negative(sleep_linear_velocity_threshold_m_s, real(0.01));
    non_negative(sleep_angular_velocity_threshold_rad_s, real(0.05));
    sleep_time_threshold_tu = std::max<tick_t>(1, sleep_time_threshold_tu);
}

} // namespace sqp
