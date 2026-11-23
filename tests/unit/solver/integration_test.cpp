#include "solver/integration.h"

#include <catch_amalgamated.hpp>

TEST_CASE("split integration matches semi-implicit body integration", "[solver][integration]") {
    sqp::Vec2 position{1, 2}, velocity{3, -4};
    sqp::real angle = sqp::real(0.25), angular_velocity = 2;
    const auto duration = sqp::config::TU_PER_SECOND / 2;

    sqp::Vec2 combined_position = position, combined_velocity = velocity;
    sqp::real combined_angle = angle, combined_angular_velocity = angular_velocity;
    sqp::integrate_body_velocity(velocity, angular_velocity, {2, 6}, 4, duration);
    CHECK(position == sqp::Vec2{1, 2});
    CHECK(angle == sqp::real(0.25));
    sqp::integrate_body_position(position, angle, velocity, angular_velocity, duration);
    sqp::integrate_body_state(combined_position, combined_angle, combined_velocity,
                              combined_angular_velocity, {2, 6}, 4, duration);

    CHECK(position == combined_position);
    CHECK(angle == combined_angle);
    CHECK(velocity == combined_velocity);
    CHECK(angular_velocity == combined_angular_velocity);
    CHECK(position == sqp::Vec2{3, sqp::real(1.5)});
    CHECK(angle == sqp::real(2.25));
}

TEST_CASE("body prediction converts simulation time units to seconds",
          "[solver][integration][prediction]") {
    sqp::Vec2 position;
    sqp::real angle = 0;
    sqp::Vec2 velocity;
    sqp::real angular_velocity = 0;

    sqp::predict_body_state({1, 2}, sqp::real(0.25), {3, -4}, 2, {}, 0,
                            sqp::config::TU_PER_SECOND / 2, position, angle, velocity,
                            angular_velocity);
    CHECK(position.x == Catch::Approx(2.5));
    CHECK(position.y == Catch::Approx(0));
    CHECK(angle == Catch::Approx(1.25));
    CHECK(velocity.x == Catch::Approx(3));
    CHECK(velocity.y == Catch::Approx(-4));
    CHECK(angular_velocity == Catch::Approx(2));

    sqp::predict_body_state({}, 0, {1, 0}, 0, {}, 0, sqp::config::TU_PER_TICK / 2, position, angle);
    CHECK(position.x == Catch::Approx(1.0 / 128.0));
    CHECK(angle == Catch::Approx(0));
}
