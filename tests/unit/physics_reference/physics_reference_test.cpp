// Supported physics reference scenarios ported from refactor's
// tests/unit/physics reference/basic_events_test.cpp to the current World API.

// Most important test file, makes sure that my physics engine "runs" normally.
// So I don't have to check every function in the visualizer manually.

#include "reference_helper.h"
#include <cmath>

TEST_CASE("Reference free motion, gravity and kinematic registration", "[physics_reference]") {
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    SECTION("stationary") {
        ReferenceWorld reference(settings);
        sqp::RigidBodyDef definition;
        definition.position = {3, 5};
        const auto body = reference.rigid(definition);
        reference.ticks(64);
        CHECK(reference.state(body).position.x == Catch::Approx(3));
        CHECK(reference.state(body).position.y == Catch::Approx(5));
    }
    SECTION("constant velocity") {
        ReferenceWorld reference(settings);
        sqp::RigidBodyDef definition;
        definition.position = {2, -1};
        definition.linear_velocity = {3, -2};
        const auto body = reference.rigid(definition);
        reference.ticks(64);
        CHECK(reference.state(body).position.x == Catch::Approx(5));
        CHECK(reference.state(body).position.y == Catch::Approx(-3));
    }
    SECTION("kinematic") {
        settings.gravity = {0, -10};
        ReferenceWorld reference(settings);
        sqp::RigidBodyDef definition;
        definition.position = {1, 3};
        definition.linear_velocity = {2, 1};
        definition.configs.register_type = sqp::BodyConfigs::RegisterType::Kinematic;
        const auto body = reference.rigid(definition);
        reference.ticks(64);
        CHECK(reference.state(body).position.x == Catch::Approx(3));
        CHECK(reference.state(body).position.y == Catch::Approx(4));
    }
    SECTION("immutable static") {
        ReferenceWorld reference(settings);
        sqp::RigidBodyDef definition;
        definition.position = {4, 1};
        definition.linear_velocity = {4, -10};
        definition.configs.register_type = sqp::BodyConfigs::RegisterType::ImmutableStatic;
        const auto body = reference.rigid(definition);
        reference.ticks(64);
        CHECK(reference.state(body).position.x == Catch::Approx(4));
        CHECK(reference.state(body).position.y == Catch::Approx(1));
        CHECK_FALSE(reference.world.query().bodies().island_of(body));
    }
}

TEST_CASE("Reference gravity, overrides and locked axes", "[physics_reference]") {
    sqp::PhysicsSettings settings;
    settings.gravity = {0, -10};
    SECTION("free fall") {
        ReferenceWorld reference(settings);
        sqp::RigidBodyDef definition;
        definition.position = {0, 10};
        const auto body = reference.rigid(definition);
        reference.ticks(64);
        CHECK(reference.state(body).linear_velocity.y == Catch::Approx(-10));
        CHECK(reference.state(body).position.y == Catch::Approx(5).margin(0.1));
    }
    SECTION("gravity override") {
        ReferenceWorld reference(settings);
        sqp::RigidBodyDef definition;
        definition.position = {0, 3};
        definition.linear_velocity = {1, 2};
        definition.configs.override_flags = sqp::BodyConfigs::OverrideFlags::NoGravity;
        const auto body = reference.rigid(definition);
        reference.ticks(64);
        CHECK(reference.state(body).position.x == Catch::Approx(1));
        CHECK(reference.state(body).position.y == Catch::Approx(5));
    }
    SECTION("locked linear axis") {
        ReferenceWorld reference(settings);
        sqp::RigidBodyDef definition;
        definition.position = {0, 3};
        definition.linear_velocity = {1, 0};
        definition.configs.override_flags = sqp::BodyConfigs::OverrideFlags::LockLinearY;
        const auto body = reference.rigid(definition);
        reference.ticks(64);
        CHECK(reference.state(body).position.x == Catch::Approx(1));
        CHECK(reference.state(body).position.y == Catch::Approx(3));
    }
    SECTION("horizontal gravity") {
        settings.gravity = {4, 0};
        ReferenceWorld reference(settings);
        const auto body = reference.rigid({});
        reference.ticks(64);
        CHECK(reference.state(body).linear_velocity.x == Catch::Approx(4));
        CHECK(reference.state(body).position.x == Catch::Approx(2).margin(0.1));
        CHECK(reference.state(body).position.y == Catch::Approx(0));
    }
}

TEST_CASE("Reference constant force and impulse", "[physics_reference]") {
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    ReferenceWorld reference(settings);
    const auto a = reference.rigid({});
    sqp::BodyProperties properties;
    REQUIRE(reference.world.query().bodies().get_properties(a, properties));
    const auto force = reference.world.set_constant_force(0, a, {properties.mass * 4, 0});
    REQUIRE(reference.world.step(0));
    reference.take(force);
    reference.ticks(64);
    CHECK(reference.state(a).linear_velocity.x == Catch::Approx(4));
    CHECK(reference.state(a).position.x == Catch::Approx(2).margin(0.1));

    const auto b = reference.rigid({});
    REQUIRE(reference.world.query().bodies().get_properties(b, properties));
    const auto impulse =
        reference.world.apply_impulse(reference.world.time_tu(), b, {properties.mass * 3, 0});
    REQUIRE(reference.world.step(0));
    reference.take(impulse);
    reference.ticks(1);
    CHECK(reference.state(b).linear_velocity.x == Catch::Approx(3));
    CHECK(reference.state(b).position.x == Catch::Approx(3.0 / 64.0).margin(1e-4));
}

TEST_CASE("Reference equal-mass head-on collisions", "[physics_reference][contact]") {
    const sqp::real restitution = GENERATE(sqp::real(0), sqp::real(0.5), sqp::real(1));
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    ReferenceWorld reference(settings);
    reference.material(restitution);
    sqp::RigidBodyDef left, right;
    left.position = {-1.5, 0};
    left.linear_velocity = {4, 0};
    right.position = {1.5, 0};
    right.linear_velocity = {-4, 0};
    const auto a = reference.rigid(left);
    const auto b = reference.rigid(right);
    reference.ticks(32);
    const auto state_a = reference.state(a);
    const auto state_b = reference.state(b);
    const double expected = 4 * restitution;
    CHECK(double(state_a.linear_velocity.x) == Catch::Approx(-expected).margin(0.05));
    CHECK(double(state_b.linear_velocity.x) == Catch::Approx(expected).margin(0.05));
    CHECK(double(state_a.linear_velocity.x + state_b.linear_velocity.x) ==
          Catch::Approx(0).margin(0.01));
}

TEST_CASE("Reference collision filtering permits uninterrupted travel",
          "[physics_reference][contact]") {
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    ReferenceWorld reference(settings);
    sqp::RigidBodyDef moving, fixed;
    moving.position = {-2, 0};
    moving.linear_velocity = {8, 0};
    fixed.configs.register_type = sqp::BodyConfigs::RegisterType::ImmutableStatic;
    reference.rigid(fixed);
    const auto body = reference.rigid(moving);
    sqp::CollisionFilter filter;
    filter.category_bits = 1;
    filter.collision_mask = 0;
    const auto command = reference.world.set_body_collision_filter(0, body, filter);
    REQUIRE(reference.world.step(0));
    reference.take(command);
    reference.ticks(32);
    CHECK(reference.state(body).position.x == Catch::Approx(2));
    CHECK(reference.state(body).linear_velocity.x == Catch::Approx(8));
}

TEST_CASE("Reference rotation, angular impulse and torque", "[physics_reference]") {
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    SECTION("uniform angular velocity") {
        ReferenceWorld reference(settings);
        sqp::RigidBodyDef definition;
        definition.linear_velocity = {1, 0};
        definition.angular_velocity = 2;
        const auto body = reference.rigid(definition);
        reference.ticks(64);
        CHECK(reference.state(body).rotation == Catch::Approx(2));
        CHECK(reference.state(body).angular_velocity == Catch::Approx(2));
    }
    SECTION("angular impulse") {
        ReferenceWorld reference(settings);
        const auto body = reference.rigid({});
        sqp::BodyProperties properties;
        REQUIRE(reference.world.query().bodies().get_properties(body, properties));
        const auto command = reference.world.apply_impulse(0, body, {}, properties.inertia * 3);
        REQUIRE(reference.world.step(0));
        reference.take(command);
        reference.ticks(1);
        CHECK(reference.state(body).angular_velocity == Catch::Approx(3));
        CHECK(reference.state(body).rotation == Catch::Approx(3.0 / 64));
    }
    SECTION("constant torque") {
        ReferenceWorld reference(settings);
        const auto body = reference.rigid({});
        sqp::BodyProperties properties;
        REQUIRE(reference.world.query().bodies().get_properties(body, properties));
        const auto command =
            reference.world.set_constant_force(0, body, {}, properties.inertia * 2);
        REQUIRE(reference.world.step(0));
        reference.take(command);
        reference.ticks(64);
        CHECK(reference.state(body).angular_velocity == Catch::Approx(2));
        CHECK(reference.state(body).rotation == Catch::Approx(1).margin(0.1));
    }
}

TEST_CASE("Reference CCD prevents a fast circle crossing a thin wall", "[physics_reference][ccd]") {
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    ReferenceWorld reference(settings);
    reference.material(1);
    sqp::RigidBodyDef wall;
    wall.configs.register_type = sqp::BodyConfigs::RegisterType::ImmutableStatic;
    reference.box(wall, {sqp::real(0.02), 10});
    sqp::RigidBodyDef moving;
    moving.position = {-5, 0};
    moving.linear_velocity = {1000, 0};
    const auto body = reference.rigid(moving, sqp::real(0.25));
    reference.ticks(1);
    CHECK(reference.state(body).position.x < sqp::real(-0.26));
    CHECK(reference.state(body).linear_velocity.x < sqp::real(-500));
}

TEST_CASE("Reference angle of repose follows the static friction threshold",
          "[physics_reference][contact][friction]") {
    constexpr sqp::real static_friction = sqp::real(0.5);
    const sqp::real repose_angle = std::atan(static_friction);
    const sqp::real degrees = sqp::math::PI / sqp::real(180);

    for (const sqp::real angle : {sqp::real(20) * degrees, sqp::real(40) * degrees}) {
        INFO("slope angle (rad): " << angle);
        sqp::PhysicsSettings settings;
        settings.gravity = {0, -10};
        ReferenceWorld reference(settings);
        reference.material(0, static_friction, sqp::real(0.3));

        const sqp::Vec2 normal{-std::sin(angle), std::cos(angle)};
        const sqp::Vec2 uphill{std::cos(angle), std::sin(angle)};
        sqp::RigidBodyDef slope;
        slope.position = normal * sqp::real(-0.5);
        slope.rotation = angle;
        slope.configs.register_type = sqp::BodyConfigs::RegisterType::ImmutableStatic;
        reference.box(slope, {10, sqp::real(0.5)});

        sqp::RigidBodyDef block;
        block.position = normal * sqp::real(0.5);
        block.rotation = angle;
        block.configs.override_flags = sqp::BodyConfigs::OverrideFlags::LockRot;
        const auto body = reference.box(block, {sqp::real(0.5), sqp::real(0.5)});
        reference.ticks(128);

        const auto state = reference.state(body);
        const sqp::real displacement =
            sqp::math::dot(state.position - block.position, uphill);
        const sqp::real speed = sqp::math::dot(state.linear_velocity, uphill);
        if (angle < repose_angle) {
            CHECK(std::abs(displacement) < sqp::real(0.1));
            CHECK(std::abs(speed) < sqp::real(0.15));
        } else {
            CHECK(displacement < sqp::real(-0.5));
            CHECK(speed < sqp::real(-0.5));
        }
    }
}

TEST_CASE("Reference three-box stack holds its horizontal alignment with one contact iteration",
          "[physics_reference][contact][stack]") {
    sqp::PhysicsSettings settings;
    settings.gravity = {0, -10};
    settings.contact.velocity_iterations = 1;
    ReferenceWorld reference(settings);
    sqp::RigidBodyDef ground;
    ground.position = {0, sqp::real(-0.5)};
    ground.configs.register_type = sqp::BodyConfigs::RegisterType::ImmutableStatic;
    reference.box(ground, {5, sqp::real(0.5)});
    std::vector<sqp::BodyHandle> stack;
    for (int row = 0; row < 3; ++row) {
        sqp::RigidBodyDef block;
        block.position = {0, sqp::real(0.5) + sqp::real(row)};
        stack.push_back(reference.box(block, {sqp::real(0.5), sqp::real(0.5)}));
    }
    reference.ticks(640);
    for (const auto body : stack)
        CHECK(std::abs(reference.state(body).position.x) < sqp::real(0.05));
}

TEST_CASE("Higher contact Hertz reduces three-box stack compression",
          "[physics_reference][contact][frequency]") {
    const auto top_sag = [](sqp::real hertz) {
        sqp::PhysicsSettings settings;
        settings.gravity = {0, -10};
        settings.contact.substeps = 4;
        settings.contact.velocity_iterations = 1;
        ReferenceWorld reference(settings);
        sqp::PhysicsMaterial material;
        material.base_restitution = 0;
        material.contact_frequency_hz = hertz;
        REQUIRE(reference.world.state().shape_pool.materials().set(0, material));
        sqp::RigidBodyDef ground;
        ground.position = {0, sqp::real(-0.5)};
        ground.configs.register_type = sqp::BodyConfigs::RegisterType::ImmutableStatic;
        reference.box(ground, {5, sqp::real(0.5)});
        sqp::BodyHandle top;
        for (int row = 0; row < 3; ++row) {
            sqp::RigidBodyDef block;
            block.position = {0, sqp::real(0.5) + sqp::real(row)};
            top = reference.box(block, {sqp::real(0.5), sqp::real(0.5)});
        }
        reference.ticks(128);
        return sqp::real(2.5) - reference.state(top).position.y;
    };
    const sqp::real soft = top_sag(16);
    const sqp::real stiff = top_sag(64);
    CHECK(stiff < soft);
}

TEST_CASE("Reference split calls agree with one long free-motion step", "[physics_reference]") {
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    ReferenceWorld one(settings), many(settings);
    sqp::RigidBodyDef definition;
    definition.linear_velocity = {5, -3};
    const auto a = one.rigid(definition);
    const auto b = many.rigid(definition);
    one.ticks(64);
    for (int tick = 0; tick < 64; ++tick)
        many.ticks(1);
    CHECK(one.state(a).position.x == Catch::Approx(many.state(b).position.x));
    CHECK(one.state(a).position.y == Catch::Approx(many.state(b).position.y));
}

TEST_CASE("Reference half-tick settings change applies gravity only afterward",
          "[physics_reference]") {
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    ReferenceWorld reference(settings);
    const auto body = reference.rigid({});
    settings.gravity = {0, -10};
    const auto command = reference.world.set_settings(sqp::config::TU_PER_TICK / 2, settings);
    reference.ticks(1);
    reference.take(command);
    CHECK(reference.state(body).linear_velocity.y == Catch::Approx(-10.0 / 128));
}

TEST_CASE("Reference lite and rigid circles exchange collision velocity",
          "[physics_reference][contact]") {
    const bool mixed = GENERATE(false, true);
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    ReferenceWorld reference(settings);
    reference.material(mixed ? sqp::real(1) : sqp::real(0));
    sqp::BodyHandle a;
    if (mixed) {
        sqp::RigidBodyDef definition;
        definition.position = {-1.5, 0};
        definition.linear_velocity = {4, 0};
        a = reference.rigid(definition);
    } else {
        a = reference.lite({-1.5, 0}, {4, 0});
    }
    const auto b = reference.lite({1.5, 0}, {-4, 0});
    reference.ticks(32);
    const double expected = mixed ? 4 : 0;
    CHECK(double(reference.state(a).linear_velocity.x) == Catch::Approx(-expected).margin(0.05));
    CHECK(double(reference.state(b).linear_velocity.x) == Catch::Approx(expected).margin(0.05));
    CHECK(reference.state(a).position.x < reference.state(b).position.x);
}

TEST_CASE("Reference overlapping circles separate", "[physics_reference][contact]") {
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    ReferenceWorld reference(settings);
    reference.material(0);
    sqp::RigidBodyDef left, right;
    left.position = {sqp::real(-0.4), 0};
    right.position = {sqp::real(0.4), 0};
    const auto a = reference.rigid(left);
    const auto b = reference.rigid(right);
    reference.ticks(128);
    CHECK(reference.state(b).position.x - reference.state(a).position.x > sqp::real(0.94));
}

TEST_CASE("Reference unequal masses conserve collision momentum", "[physics_reference][contact]") {
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    ReferenceWorld reference(settings);
    reference.material(0);
    sqp::RigidBodyDef moving, waiting;
    moving.position = {-3, 0};
    moving.linear_velocity = {6, 0};
    const auto a = reference.rigid(moving, sqp::real(0.5));
    const auto b = reference.rigid(waiting, 1);
    sqp::BodyProperties props_a, props_b;
    REQUIRE(reference.world.query().bodies().get_properties(a, props_a));
    REQUIRE(reference.world.query().bodies().get_properties(b, props_b));
    reference.ticks(40);
    const double expected = double(props_a.mass / (props_a.mass + props_b.mass) * 6);
    CHECK(double(reference.state(a).linear_velocity.x) == Catch::Approx(expected).margin(0.05));
    CHECK(double(reference.state(b).linear_velocity.x) == Catch::Approx(expected).margin(0.05));
    CHECK(double(props_a.mass * reference.state(a).linear_velocity.x +
                 props_b.mass * reference.state(b).linear_velocity.x) ==
          Catch::Approx(double(props_a.mass * 6)).margin(0.05));
}

TEST_CASE("Reference falling circle rests on a static floor", "[physics_reference][contact]") {
    sqp::PhysicsSettings settings;
    settings.gravity = {0, -10};
    ReferenceWorld reference(settings);
    reference.material(0, sqp::real(0.8), sqp::real(0.6));
    sqp::RigidBodyDef floor, falling;
    floor.position = {0, sqp::real(-0.5)};
    floor.configs.register_type = sqp::BodyConfigs::RegisterType::ImmutableStatic;
    reference.box(floor, {15, sqp::real(0.5)});
    falling.position = {0, 2};
    const auto body = reference.rigid(falling);
    reference.ticks(128);
    CHECK(double(reference.state(body).position.y) == Catch::Approx(0.5).margin(0.02));
    CHECK(double(reference.state(body).linear_velocity.y) == Catch::Approx(0).margin(0.05));
}

TEST_CASE("Reference elastic floor bounce", "[physics_reference][contact]") {
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    ReferenceWorld reference(settings);
    reference.material(1);
    sqp::RigidBodyDef floor, falling;
    floor.position = {0, sqp::real(-0.5)};
    floor.configs.register_type = sqp::BodyConfigs::RegisterType::ImmutableStatic;
    reference.box(floor, {10, sqp::real(0.5)});
    falling.position = {0, 1};
    falling.linear_velocity = {0, -4};
    const auto body = reference.rigid(falling);
    reference.ticks(32);
    CHECK(double(reference.state(body).linear_velocity.y) == Catch::Approx(4).margin(0.05));
    CHECK(double(reference.state(body).position.y) == Catch::Approx(2).margin(0.1));
}

TEST_CASE("Reference oblique impact preserves common tangential velocity",
          "[physics_reference][contact]") {
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    ReferenceWorld reference(settings);
    reference.material(1);
    sqp::RigidBodyDef left, right;
    left.position = {-1.5, 0};
    left.linear_velocity = {4, 1};
    right.position = {1.5, 0};
    right.linear_velocity = {-4, 1};
    const auto a = reference.rigid(left);
    const auto b = reference.rigid(right);
    reference.ticks(32);
    CHECK(double(reference.state(a).linear_velocity.x) == Catch::Approx(-4).margin(0.05));
    CHECK(double(reference.state(b).linear_velocity.x) == Catch::Approx(4).margin(0.05));
    CHECK(double(reference.state(a).linear_velocity.y) == Catch::Approx(1).margin(0.05));
    CHECK(double(reference.state(b).linear_velocity.y) == Catch::Approx(1).margin(0.05));
}

TEST_CASE("Reference fast lite circle reflects from a thin wall", "[physics_reference][ccd]") {
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    ReferenceWorld reference(settings);
    reference.material(1);
    sqp::RigidBodyDef wall;
    wall.configs.register_type = sqp::BodyConfigs::RegisterType::ImmutableStatic;
    reference.box(wall, {sqp::real(0.02), 10});
    const auto body = reference.lite({-5, 0}, {1000, 0}, sqp::real(0.25));
    reference.ticks(1);
    CHECK(reference.state(body).position.x < sqp::real(-0.26));
    CHECK(reference.state(body).linear_velocity.x < sqp::real(-500));
}

TEST_CASE("Reference collision conserves momentum and follows restitution energy",
          "[physics_reference][contact]") {
    const sqp::real restitution = GENERATE(sqp::real(0), sqp::real(1));
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    ReferenceWorld reference(settings);
    reference.material(restitution);
    sqp::RigidBodyDef left, right;
    left.position = {-1.5, 0};
    left.linear_velocity = {4, 0};
    right.position = {1.5, 0};
    right.linear_velocity = {-4, 0};
    const auto a = reference.rigid(left);
    const auto b = reference.rigid(right);
    sqp::BodyProperties pa, pb;
    REQUIRE(reference.world.query().bodies().get_properties(a, pa));
    REQUIRE(reference.world.query().bodies().get_properties(b, pb));
    const double initial_energy = double((pa.mass + pb.mass) * 8);
    reference.ticks(32);
    const auto va = reference.state(a).linear_velocity;
    const auto vb = reference.state(b).linear_velocity;
    const double momentum = double(pa.mass * va.x + pb.mass * vb.x);
    const double final_energy =
        double((pa.mass * va.length_squared() + pb.mass * vb.length_squared()) / 2);
    CHECK(momentum == Catch::Approx(0).margin(0.05));
    if (restitution == sqp::real(0))
        CHECK(final_energy < initial_energy * 0.05);
    else
        CHECK(final_energy == Catch::Approx(initial_energy).margin(0.05));
}

TEST_CASE("Reference floor contact preserves tangential travel with CCD modes",
          "[physics_reference][ccd][contact]") {
    const auto mode =
        GENERATE(sqp::BodyConfigs::CcdOverride::Auto, sqp::BodyConfigs::CcdOverride::Disabled,
                 sqp::BodyConfigs::CcdOverride::Shapecast);
    sqp::PhysicsSettings settings;
    settings.gravity = {0, -10};
    ReferenceWorld reference(settings);
    reference.material(0);
    sqp::RigidBodyDef floor, moving;
    floor.position = {0, sqp::real(-0.5)};
    floor.configs.register_type = sqp::BodyConfigs::RegisterType::ImmutableStatic;
    reference.box(floor, {20, sqp::real(0.5)});
    moving.position = {0, sqp::real(0.5)};
    moving.linear_velocity = {2, 0};
    const auto body = reference.rigid(moving);
    const auto command = reference.world.set_body_ccd_override(0, body, mode);
    REQUIRE(reference.world.step(0));
    reference.take(command);
    reference.ticks(64);
    CHECK(double(reference.state(body).position.x) == Catch::Approx(2).margin(0.05));
    CHECK(double(reference.state(body).linear_velocity.x) == Catch::Approx(2).margin(0.05));
}

TEST_CASE("Reference rolling circle continues along the floor",
          "[physics_reference][ccd][contact]") {
    const auto mode =
        GENERATE(sqp::BodyConfigs::CcdOverride::Auto, sqp::BodyConfigs::CcdOverride::Disabled);
    sqp::PhysicsSettings settings;
    settings.gravity = {0, -10};
    ReferenceWorld reference(settings);
    reference.material(0, sqp::real(0.8), sqp::real(0.8));
    sqp::RigidBodyDef floor, moving;
    floor.position = {0, sqp::real(-0.5)};
    floor.configs.register_type = sqp::BodyConfigs::RegisterType::ImmutableStatic;
    reference.box(floor, {20, sqp::real(0.5)});
    moving.position = {0, sqp::real(0.5)};
    moving.linear_velocity = {4, 0};
    const auto body = reference.rigid(moving);
    const auto command = reference.world.set_body_ccd_override(0, body, mode);
    REQUIRE(reference.world.step(0));
    reference.take(command);
    reference.ticks(64);
    CHECK(reference.state(body).position.x > sqp::real(2));
}
