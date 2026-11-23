#include "sqphysics2d/settings.h"

#include <catch_amalgamated.hpp>

#include <limits>

TEST_CASE("interval settings normalize one minimum period and TOI budget", "[world][settings]") {
    sqp::IntervalSettings settings;
    settings.min_period_tu = 0;
    settings.max_truncations_per_tick = 0;
    settings.sanitize();
    CHECK(settings.min_period_tu == sqp::config::TU_PER_TICK / 4);
    CHECK(settings.max_truncations_per_tick == 1);

    settings.min_period_tu = 3;
    settings.sanitize();
    CHECK(settings.min_period_tu == 4);

    settings.min_period_tu = sqp::config::TU_PER_TICK + 1;
    settings.sanitize();
    CHECK(settings.min_period_tu == sqp::config::TU_PER_TICK);
}

TEST_CASE("physics settings sanitize invalid solver and motion values", "[world][settings]") {
    sqp::PhysicsSettings settings;
    settings.interval.max_linear_displacement_m = -1;
    settings.interval.max_angular_displacement_rad = std::numeric_limits<sqp::real>::infinity();
    settings.contact.substeps = 0;
    settings.contact.velocity_iterations = 300;
    settings.contact.toi_velocity_iterations = 0;
    settings.contact.restitution_threshold_m_s = -1;
    settings.ccd.linear_start_threshold_m = -1;
    settings.broad_phase.proxy_extension_ratio = std::numeric_limits<sqp::real>::infinity();
    settings.gravity.x = std::numeric_limits<sqp::real>::infinity();
    settings.linear_velocity_tick_multiplier = -1;
    settings.angular_velocity_tick_multiplier = std::numeric_limits<sqp::real>::quiet_NaN();
    settings.sleep_linear_velocity_threshold_m_s = -1;
    settings.sleep_angular_velocity_threshold_rad_s = std::numeric_limits<sqp::real>::infinity();
    settings.sleep_time_threshold_tu = 0;

    settings.sanitize();

    CHECK(settings.interval.max_linear_displacement_m == sqp::real(0.1));
    CHECK(settings.interval.max_angular_displacement_rad == sqp::real(0.01));
    CHECK(settings.contact.substeps == 1);
    CHECK(settings.contact.velocity_iterations == 256);
    CHECK(settings.contact.toi_velocity_iterations == 1);
    CHECK(settings.contact.restitution_threshold_m_s == sqp::real(0.1));
    CHECK(settings.ccd.linear_start_threshold_m == sqp::real(0.4));
    CHECK(settings.broad_phase.proxy_extension_ratio ==
          sqp::BroadPhaseSettings::default_proxy_extension_ratio);
    CHECK(settings.gravity.x == 0);
    CHECK(settings.gravity.y == sqp::real(-9.81));
    CHECK(settings.linear_velocity_tick_multiplier == 1);
    CHECK(settings.angular_velocity_tick_multiplier == 1);
    CHECK(settings.sleep_linear_velocity_threshold_m_s == sqp::real(0.01));
    CHECK(settings.sleep_angular_velocity_threshold_rad_s == sqp::real(0.05));
    CHECK(settings.sleep_time_threshold_tu == 1);
}

TEST_CASE("broad phase ratio clamps to its supported range", "[world][settings]") {
    sqp::PhysicsSettings settings;
    settings.broad_phase.proxy_extension_ratio = -1;
    settings.sanitize();
    CHECK(settings.broad_phase.proxy_extension_ratio == 0);
    settings.broad_phase.proxy_extension_ratio = 2;
    settings.sanitize();
    CHECK(settings.broad_phase.proxy_extension_ratio == 1);
}
