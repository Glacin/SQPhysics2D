#include "reference_helper.h"
#include <cmath>

// Tests that check friction is applied correctly.


#include "reference_helper.h"
#include <cmath>

namespace {

sqp::real simulate_friction(sqp::real friction, int ticks = 8,
                            sqp::real velocity = 8,
                            sqp::real gravity = -10,
                            sqp::real height = sqp::real(0.5)) {
    sqp::PhysicsSettings settings;
    settings.gravity = {0, gravity};
    ReferenceWorld reference(settings);
    reference.material(0, friction, friction);

    sqp::RigidBodyDef floor;
    floor.position = {0, sqp::real(-0.5)};
    floor.configs.register_type = sqp::BodyConfigs::RegisterType::ImmutableStatic;
    reference.box(floor, {20, sqp::real(0.5)});

    sqp::RigidBodyDef moving;
    moving.position = {0, height};
    moving.linear_velocity = {velocity, 0};
    const auto body = reference.rigid(moving);

    reference.ticks(ticks);
    return reference.state(body).linear_velocity.x;
}

} // namespace

TEST_CASE("Reference contact friction follows Coulomb sliding law",
          "[physics_reference][contact][friction]") {
    // v(t) = v0 - mu * g * t
    // t = 8 / 64 = 0.125s
    const auto v0 = simulate_friction(0);
    const auto v025 = simulate_friction(sqp::real(0.25));
    const auto v050 = simulate_friction(sqp::real(0.5));
    const auto v100 = simulate_friction(1);

    CHECK(std::abs(v0 - sqp::real(8)) < sqp::real(0.03));
    CHECK(std::abs(v025 - sqp::real(7.6875)) < sqp::real(0.2));
    CHECK(std::abs(v050 - sqp::real(7.375)) < sqp::real(0.2));
    CHECK(std::abs(v100 - sqp::real(6.75)) < sqp::real(0.2));

    CHECK(v0 > v025);
    CHECK(v025 > v050);
    CHECK(v050 > v100);
}

TEST_CASE("Reference contact friction is symmetric in sliding direction",
          "[physics_reference][contact][friction]") {
    const auto positive = simulate_friction(sqp::real(0.5), 8, 8);
    const auto negative = simulate_friction(sqp::real(0.5), 8, -8);

    CHECK(positive > 0);
    CHECK(negative < 0);
    CHECK(std::abs(positive + negative) < sqp::real(0.03));
}

TEST_CASE("Reference contact friction scales with normal force",
          "[physics_reference][contact][friction]") {
    const auto weak = simulate_friction(sqp::real(0.5), 8, 8, -5);
    const auto normal = simulate_friction(sqp::real(0.5), 8, 8, -10);
    const auto strong = simulate_friction(sqp::real(0.5), 8, 8, -20);

    CHECK(weak > normal);
    CHECK(normal > strong);

    CHECK(std::abs(weak - sqp::real(7.6875)) < sqp::real(0.2));
    CHECK(std::abs(normal - sqp::real(7.375)) < sqp::real(0.2));
    CHECK(std::abs(strong - sqp::real(6.75)) < sqp::real(0.2));
}

TEST_CASE("Reference contact friction requires contact",
          "[physics_reference][contact][friction]") {
    // The body remains well above the floor throughout the simulation.
    const auto free = simulate_friction(0, 8, 8, -10, 5);
    const auto friction = simulate_friction(1, 8, 8, -10, 5);

    CHECK(std::abs(free - sqp::real(8)) < sqp::real(0.03));
    CHECK(std::abs(friction - sqp::real(8)) < sqp::real(0.03));
    CHECK(std::abs(free - friction) < sqp::real(0.01));
}

TEST_CASE("Reference contact friction does not reverse sliding direction",
          "[physics_reference][contact][friction]") {
    for (const auto friction : {sqp::real(0.25), sqp::real(0.5), sqp::real(1)}) {
        const auto positive = simulate_friction(friction, 8, 8);
        const auto negative = simulate_friction(friction, 8, -8);

        CHECK(positive >= 0);
        CHECK(positive < 8);
        CHECK(negative <= 0);
        CHECK(negative > -8);
    }
}

TEST_CASE("Reference resting contact does not generate tangential velocity",
          "[physics_reference][contact][friction]") {
    for (const auto friction : {sqp::real(0), sqp::real(0.5), sqp::real(1)}) {
        const auto velocity = simulate_friction(friction, 128, 0);

        CHECK(std::abs(velocity) < sqp::real(0.01));
    }
}
