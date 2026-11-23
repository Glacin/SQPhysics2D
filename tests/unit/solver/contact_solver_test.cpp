#include <algorithm>
#include <cmath>
#include <limits>

#include "solver/contact_material.h"

#include "sqphysics2d/detail/contact_solver.h"
#include "sqphysics2d/detail/data_storage.h"
#include <catch_amalgamated.hpp>

namespace {
sqp::ContactConstraint prepare_from_materials(const sqp::BodyPool& bodies, sqp::BodyHandle a,
                                              sqp::BodyHandle b,
                                              const sqp::ContactManifold& manifold,
                                              const sqp::PhysicsMaterial& first,
                                              const sqp::PhysicsMaterial& second, sqp::real scale_a,
                                              sqp::real scale_b, sqp::real h, bool persisted) {
    const auto combined = sqp::detail::combine_contact_materials(first, second);
    sqp::ContactSolverParameters parameters{
        combined.restitution,          combined.dynamic_friction,
        combined.static_friction,      combined.contact_damping_ratio,
        combined.contact_frequency_hz, combined.max_penetration_correction_speed};
    sqp::BodyState state_a, state_b;
    sqp::BodyProperties properties_a, properties_b;
    if (manifold.count > 0 && bodies.get_state(a, state_a) && bodies.get_state(b, state_b) &&
        bodies.get_properties(a, properties_a) && bodies.get_properties(b, properties_b)) {
        const auto length = manifold.normal.length();
        if (length > 0) {
            const auto normal = manifold.normal / length;
            const auto point = manifold.points[0].point;
            const auto closing = -sqp::math::dot(state_b.velocity_at_world_point(point) -
                                                     state_a.velocity_at_world_point(point),
                                                 normal);
            parameters.restitution = sqp::detail::contact_restitution(
                combined, closing, properties_a.inv_mass + properties_b.inv_mass,
                std::min(scale_a, scale_b), true);
        }
    }
    return sqp::prepare_contact_constraint(bodies, a, b, manifold, parameters, h, persisted);
}

sqp::BodyHandle
add_rigid(sqp::BodyPool& bodies, sqp::ShapePool& shapes, sqp::Vec2 position, sqp::Vec2 velocity,
          sqp::BodyConfigs::RegisterType type = sqp::BodyConfigs::RegisterType::Dynamic,
          sqp::u32 overrides = 0) {
    sqp::RigidBodyDef def(position, 0, velocity, 0);
    def.configs.register_type = type;
    def.configs.override_flags = sqp::BodyConfigs::OverrideFlags(overrides);
    const sqp::BodyHandle body = bodies.add_rigid(def);
    const sqp::ShapeHandle shape =
        shapes.add_circle({sqp::real(1), {}, 0, sqp::real(1) / sqp::math::PI});
    REQUIRE(bodies.attach_shape(body, shape, shapes));
    return body;
}

sqp::BodyHandle add_lite(sqp::BodyPool& bodies, sqp::ShapePool& shapes, sqp::Vec2 position,
                         sqp::Vec2 velocity, sqp::u32 overrides = 0) {
    const sqp::ShapeHandle shape =
        shapes.add_circle({sqp::real(1), {}, 0, sqp::real(1) / sqp::math::PI});
    sqp::LiteCircleBodyDef def(position, velocity);
    def.configs.override_flags = sqp::BodyConfigs::OverrideFlags(overrides);
    return bodies.add_lite_circle(def, shape, shapes);
}

sqp::ContactManifold contact(sqp::Vec2 point, sqp::Vec2 normal, sqp::real depth) {
    sqp::ContactManifold manifold;
    manifold.count = 1;
    manifold.normal = normal;
    manifold.points[0].point = point;
    manifold.points[0].depth = depth;
    return manifold;
}

sqp::Vec2 velocity(const sqp::BodyPool& bodies, sqp::BodyHandle body) {
    sqp::BodyState state;
    REQUIRE(bodies.get_state(body, state));
    return state.linear_velocity;
}
} // namespace

TEST_CASE("Contact softness handles disabled and finite frequencies", "[contact][softness]") {
    CHECK(sqp::make_contact_softness(0, 1, sqp::real(0.01)).mass_scale == 0);
    CHECK(sqp::make_contact_softness(10, 1, 0).bias_rate == 0);
    const auto softness = sqp::make_contact_softness(10, 2, sqp::real(0.01));
    CHECK(std::isfinite(softness.bias_rate));
    CHECK(softness.bias_rate > 0);
    CHECK(softness.mass_scale > 0);
    CHECK(softness.impulse_scale > 0);
    const auto stiff = sqp::make_contact_softness(std::numeric_limits<sqp::real>::max(),
                                                 8, sqp::real(1) / 256);
    CHECK(std::isfinite(stiff.bias_rate));
    CHECK(stiff.mass_scale == 1);
    CHECK(stiff.impulse_scale == 0);
}

TEST_CASE("Contact preparation uses all 64 material Hertz in four substeps",
          "[contact][softness][frequency]") {
    sqp::BodyPool bodies;
    sqp::ShapePool shapes;
    const auto first = add_rigid(bodies, shapes, {}, {});
    const auto second = add_rigid(bodies, shapes, {0, 1}, {});
    sqp::ContactSolverParameters parameters;
    parameters.contact_frequency = 64;
    parameters.damping_ratio = 8;
    const auto manifold = contact({0, sqp::real(0.5)}, {0, 1}, sqp::real(0.01));
    const sqp::real h = sqp::real(1) / 256;
    const auto dynamic = sqp::prepare_contact_constraint(bodies, first, second, manifold,
                                                          parameters, h, false);
    const auto expected_dynamic = sqp::make_contact_softness(64, 8, h);
    CHECK(dynamic.softness.mass_scale == Catch::Approx(expected_dynamic.mass_scale));
    CHECK(dynamic.softness.mass_scale > sqp::make_contact_softness(32, 8, h).mass_scale);

    const auto ground = add_rigid(bodies, shapes, {}, {},
                                  sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto static_contact = sqp::prepare_contact_constraint(bodies, ground, second,
                                                                 manifold, parameters, h, false);
    const auto expected_static = sqp::make_contact_softness(128, 8, h);
    CHECK(static_contact.softness.mass_scale == Catch::Approx(expected_static.mass_scale));
}

TEST_CASE("Box2D normal PGS vectors agree for one and two points",
          "[contact][pgs][box2d]") {
    struct Case {
        bool static_a;
        int point_count;
        sqp::real depth;
        sqp::real closing_speed;
        sqp::real expected_impulse;
        sqp::real expected_velocity_b;
    };
    // Reference values follow Box2D's sequential normal recurrence with fixed
    // contact arms, one push and one relax pass, and rigid mass scaling.
    const Case cases[] = {{true, 1, 0, 1, 1, 0},
                          {true, 2, 0, 1, 1, 0},
                          {false, 2, 0, 1, sqp::real(0.5), sqp::real(-0.5)},
                          {true, 1, sqp::real(-0.01), 2, 1, -1}};
    for (const auto& example : cases) {
        sqp::BodyPool bodies;
        sqp::ShapePool shapes;
        const auto flags = static_cast<sqp::u32>(sqp::BodyConfigs::OverrideFlags::LockRot);
        const auto a = add_rigid(bodies, shapes, {}, {},
                                 example.static_a
                                     ? sqp::BodyConfigs::RegisterType::ImmutableStatic
                                     : sqp::BodyConfigs::RegisterType::Dynamic,
                                 flags);
        const auto b = add_rigid(bodies, shapes, {0, 1}, {0, -example.closing_speed},
                                 sqp::BodyConfigs::RegisterType::Dynamic, flags);
        auto manifold = contact({sqp::real(-0.5), sqp::real(0.5)}, {0, 1}, example.depth);
        manifold.count = example.point_count;
        manifold.points[1] = manifold.points[0];
        manifold.points[1].point.x = sqp::real(0.5);
        sqp::ContactSolverParameters parameters;
        auto constraint = sqp::prepare_contact_constraint(bodies, a, b, manifold,
                                                           parameters, sqp::real(0.01), true);
        constraint.softness = {0, 1, 0};
        sqp::PreparedContactConstraint prepared;
        REQUIRE(sqp::prepare_contact_runtime(bodies, constraint, prepared));
        sqp::solve_contact_constraint(prepared, 0, true);
        sqp::solve_contact_constraint(prepared, 0, false);
        CHECK(constraint.points[0].normal_impulse ==
              Catch::Approx(example.expected_impulse));
        if (example.point_count == 2)
            CHECK(constraint.points[1].normal_impulse == Catch::Approx(0));
        CHECK(velocity(bodies, b).y == Catch::Approx(example.expected_velocity_b));
    }
}

TEST_CASE("Box2D two-point normal recurrence includes angular cross coupling",
          "[contact][pgs][box2d]") {
    sqp::BodyPool bodies;
    sqp::ShapePool shapes;
    const auto ground = add_rigid(bodies, shapes, {}, {},
                                  sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto body = add_rigid(bodies, shapes, {0, 1}, {0, -1});
    auto manifold = contact({sqp::real(-0.5), sqp::real(0.5)}, {0, 1}, 0);
    manifold.count = 2;
    manifold.points[1] = manifold.points[0];
    manifold.points[1].point.x = sqp::real(0.5);
    sqp::ContactSolverParameters parameters;
    auto constraint = sqp::prepare_contact_constraint(bodies, ground, body, manifold,
                                                       parameters, sqp::real(0.01), true);
    constraint.softness = {0, 1, 0};
    sqp::PreparedContactConstraint prepared;
    REQUIRE(sqp::prepare_contact_runtime(bodies, constraint, prepared));
    sqp::solve_contact_constraint(prepared, 0, true);
    CHECK(constraint.points[0].normal_impulse == Catch::Approx(sqp::real(2) / 3));
    CHECK(constraint.points[1].normal_impulse == Catch::Approx(sqp::real(4) / 9));
    sqp::solve_contact_constraint(prepared, 0, false);
    CHECK(constraint.points[0].normal_impulse == Catch::Approx(sqp::real(14) / 27));
    CHECK(constraint.points[1].normal_impulse == Catch::Approx(sqp::real(40) / 81));
    CHECK(velocity(bodies, body).y == Catch::Approx(sqp::real(1) / 81));
    sqp::BodyState state;
    REQUIRE(bodies.get_state(body, state));
    CHECK(state.angular_velocity == Catch::Approx(sqp::real(-2) / 81));
}

TEST_CASE("Contact preparation consumes supplied PreSolve parameters", "[contact][prepared]") {
    sqp::BodyPool bodies;
    sqp::ShapePool shapes;
    const auto ground =
        add_rigid(bodies, shapes, {}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto moving = add_rigid(bodies, shapes, {0, 1}, {0, -2});
    sqp::ContactSolverParameters values;
    values.restitution = sqp::real(0.35);
    values.dynamic_friction = sqp::real(0.2);
    values.static_friction = sqp::real(0.4);
    values.damping_ratio = sqp::real(3);
    values.contact_frequency = sqp::real(12);
    values.max_penetration_correction_speed = sqp::real(2);
    const auto constraint = sqp::prepare_contact_constraint(bodies, ground, moving,
                                                            contact({0, sqp::real(0.5)}, {0, 1}, 0),
                                                            values, sqp::real(0.01), false);
    CHECK(constraint.restitution == values.restitution);
    CHECK(constraint.dynamic_friction == values.dynamic_friction);
    CHECK(constraint.static_friction == values.static_friction);
    CHECK(constraint.max_penetration_correction_speed == values.max_penetration_correction_speed);
    CHECK(constraint.points[0].restitution_velocity == Catch::Approx(0.7));
    CHECK(constraint.softness.mass_scale > 0);
}

TEST_CASE("Pair-local TOI constraints retain anchors and restitution near x=1024",
          "[contact][toi]") {
    for (const bool reverse : {false, true}) {
        sqp::BodyPool bodies;
        sqp::ShapePool shapes;
        const auto wall = add_rigid(bodies, shapes, {sqp::real(1024.5), 0}, {},
                                    sqp::BodyConfigs::RegisterType::ImmutableStatic);
        const auto moving = add_rigid(bodies, shapes, {sqp::real(1022.5), 0}, {10, 0});
        const auto a = reverse ? moving : wall;
        const auto b = reverse ? wall : moving;
        const sqp::Vec2 origin = reverse ? sqp::Vec2{sqp::real(1022.5), 0}
                                         : sqp::Vec2{sqp::real(1024.5), 0};
        auto manifold = contact({reverse ? sqp::real(1) : sqp::real(-1), 0},
                                {reverse ? sqp::real(1) : sqp::real(-1), 0}, 0);
        sqp::ContactSolverParameters parameters;
        parameters.restitution = 1;
        auto constraint = sqp::prepare_toi_local_contact_constraint(
            bodies, a, b, manifold, origin, parameters, sqp::real(0.01), false);
        REQUIRE(constraint.count == 1);
        CHECK(constraint.local_geometry);
        CHECK(constraint.points[0].anchor_a.x == Catch::Approx(reverse ? 1 : -1));
        CHECK(constraint.points[0].anchor_b.x == Catch::Approx(reverse ? -1 : 1));
        sqp::PreparedContactConstraint prepared;
        REQUIRE(sqp::prepare_contact_runtime(bodies, constraint, prepared));
        sqp::solve_contact_constraint(prepared, 0, false);
        sqp::refresh_contact_constraint_geometry(prepared);
        sqp::apply_contact_restitution(prepared);
        CHECK(velocity(bodies, moving).x < 0);
        const auto separation = constraint.points[0].current_separation;
        sqp::BodyState moved;
        REQUIRE(bodies.get_state(moving, moved));
        moved.position.x += sqp::real(0.125);
        REQUIRE(bodies.set_state(moving, moved));
        sqp::refresh_contact_constraint_geometry(prepared);
        CHECK(constraint.points[0].current_separation != separation);
    }
}

TEST_CASE("Contact constraint geometry follows the body pose without changing its pose",
          "[contact][geometry]") {
    sqp::BodyPool bodies;
    sqp::ShapePool shapes;
    const auto ground =
        add_rigid(bodies, shapes, {}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto moving = add_rigid(bodies, shapes, {0, 1}, {});
    sqp::PhysicsMaterial material;
    auto constraint = prepare_from_materials(
        bodies, ground, moving, contact({sqp::real(0.4), sqp::real(0.5)}, {0, 1}, sqp::real(0.1)),
        material, material, 1, 1, sqp::real(0.01), false);
    const auto old_mass = constraint.points[0].normal_mass;
    const auto old_separation = constraint.points[0].base_separation;
    sqp::BodyState moved;
    REQUIRE(bodies.get_state(moving, moved));
    moved.position.y -= sqp::real(0.2);
    moved.rotation = sqp::real(0.5);
    REQUIRE(bodies.set_state(moving, moved));
    sqp::PreparedContactConstraint prepared;
    REQUIRE(sqp::prepare_contact_runtime(bodies, constraint, prepared));
    sqp::refresh_contact_constraint_geometry(prepared);
    CHECK(constraint.points[0].base_separation == Catch::Approx(old_separation));
    CHECK(constraint.points[0].current_separation != Catch::Approx(old_separation));
    CHECK(constraint.points[0].normal_mass == Catch::Approx(old_mass));
    sqp::BodyState after;
    REQUIRE(bodies.get_state(moving, after));
    CHECK(after == moved);
}

TEST_CASE("Prepared two-point geometry refreshes translation and rotation across substeps",
          "[contact][geometry][prepared]") {
    sqp::BodyPool bodies;
    sqp::ShapePool shapes;
    const auto ground =
        add_rigid(bodies, shapes, {}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto moving = add_rigid(bodies, shapes, {0, 1}, {});
    sqp::ContactManifold manifold =
        contact({sqp::real(-0.4), sqp::real(0.5)}, {0, 1}, sqp::real(0.1));
    manifold.count = 2;
    manifold.points[1] = manifold.points[0];
    manifold.points[1].point.x = sqp::real(0.4);
    sqp::PhysicsMaterial material;
    auto constraint = prepare_from_materials(bodies, ground, moving, manifold, material, material,
                                             1, 1, sqp::real(0.01), false);
    sqp::PreparedContactConstraint prepared;
    REQUIRE(sqp::prepare_contact_runtime(bodies, constraint, prepared));
    const sqp::real initial_mass = constraint.points[0].normal_mass;
    const sqp::real initial_anchor_y = constraint.points[0].anchor_b.y;

    sqp::BodyState moved;
    REQUIRE(bodies.get_state(moving, moved));
    moved.position.y -= sqp::real(0.2);
    REQUIRE(bodies.set_state(moving, moved));
    sqp::refresh_contact_constraint_geometry(prepared);
    for (const auto& point : constraint.points) {
        CHECK(point.current_separation == Catch::Approx(-sqp::real(0.3)));
        CHECK(point.normal_mass == initial_mass);
        CHECK(point.anchor_b.y == initial_anchor_y);
    }
    sqp::refresh_contact_constraint_geometry(prepared);
    CHECK(constraint.points[0].current_separation == Catch::Approx(-sqp::real(0.3)));

    for (sqp::real angle : {sqp::real(0.5), sqp::real(0.8)}) {
        moved.rotation = angle;
        moved.position.y -= sqp::real(0.1);
        REQUIRE(bodies.set_state(moving, moved));
        sqp::refresh_contact_constraint_geometry(prepared);
        for (sqp::i32 i = 0; i < constraint.count; ++i) {
            const sqp::real original_x = manifold.points[i].point.x;
            const sqp::real rotated_y =
                std::sin(angle) * original_x + std::cos(angle) * initial_anchor_y;
            CHECK(constraint.points[i].current_separation ==
                  Catch::Approx(-sqp::real(0.1) + moved.position.y - sqp::real(1) + rotated_y -
                                initial_anchor_y)
                      .margin(sqp::real(1e-5)));
            CHECK(constraint.points[i].anchor_b.y == Catch::Approx(initial_anchor_y));
        }
        const sqp::real separation = constraint.points[0].current_separation;
        sqp::refresh_contact_constraint_geometry(prepared);
        CHECK(constraint.points[0].current_separation == separation);
        sqp::solve_contact_constraint(prepared, 0, false);
    }
    CHECK(constraint.points[0].normal_mass == Catch::Approx(initial_mass));
}

TEST_CASE("Prepared contact endpoints are transient and restitution requires a collision impulse",
          "[contact][prepared]") {
    sqp::BodyPool bodies;
    sqp::ShapePool shapes;
    const auto ground =
        add_rigid(bodies, shapes, {0, 0}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto body = add_rigid(bodies, shapes, {0, 1}, {0, -2});
    sqp::PhysicsMaterial material;
    material.restitution_flags = sqp::RESTITUTION_DISABLE_IMPACT_SCALING;
    material.base_restitution = sqp::real(1);
    auto constraint = prepare_from_materials(bodies, ground, body,
                                             contact({0, sqp::real(0.5)}, {0, 1}, sqp::real(0.1)),
                                             material, material, 1, 1, sqp::real(0.01), false);
    sqp::PreparedContactConstraint prepared;
    REQUIRE(sqp::prepare_contact_runtime(bodies, constraint, prepared));
    sqp::apply_contact_restitution(prepared);
    CHECK(velocity(bodies, body).y == Catch::Approx(-2));
    CHECK_FALSE(constraint.restitution_applied);
    sqp::solve_contact_constraint(prepared, 0, false);
    sqp::apply_contact_restitution(prepared);
    CHECK(velocity(bodies, body).y == Catch::Approx(2));
    CHECK(constraint.restitution_applied);
}

TEST_CASE("Contact relaxation applies restitution once and reports impulses",
          "[contact][restitution]") {
    sqp::BodyPool bodies;
    sqp::ShapePool shapes;
    const auto ground =
        add_rigid(bodies, shapes, {0, 0}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto body = add_rigid(bodies, shapes, {0, 1}, {0, -2});
    sqp::PhysicsMaterial material;
    material.restitution_flags = sqp::RESTITUTION_DISABLE_IMPACT_SCALING;
    material.base_restitution = sqp::real(1);
    const auto manifold = contact({0, sqp::real(0.5)}, {0, 1}, sqp::real(0.1));
    auto constraint = prepare_from_materials(bodies, ground, body, manifold, material, material, 1,
                                             1, sqp::real(0.01), false);
    REQUIRE(constraint.count == 1);
    CHECK(constraint.points[0].restitution_velocity == Catch::Approx(2));
    sqp::PreparedContactConstraint prepared;
    REQUIRE(sqp::prepare_contact_runtime(bodies, constraint, prepared));
    sqp::solve_contact_constraint(prepared, 0, false);
    sqp::apply_contact_restitution(prepared);
    CHECK(velocity(bodies, body).y == Catch::Approx(2));
    const sqp::real impulse = constraint.points[0].normal_impulse;
    REQUIRE(impulse > 0);
    sqp::apply_contact_restitution(prepared);
    CHECK(constraint.points[0].normal_impulse == Catch::Approx(impulse));
    const auto report = sqp::report_contact_constraint(constraint);
    CHECK(report.restitution_applied);
    CHECK(report.normal_impulse > 0);
    CHECK(report.total_impulse >= report.normal_impulse);
    CHECK(constraint.normal.y == 1);
}

TEST_CASE("Contact bias removes overlap without injecting downward energy", "[contact][bias]") {
    sqp::BodyPool bodies;
    sqp::ShapePool shapes;
    const auto ground =
        add_rigid(bodies, shapes, {0, 0}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto body = add_rigid(bodies, shapes, {0, 1}, {});
    sqp::PhysicsMaterial material;
    material.contact_frequency_hz = sqp::real(20);
    auto constraint = prepare_from_materials(bodies, ground, body,
                                             contact({0, sqp::real(0.5)}, {0, 1}, sqp::real(0.2)),
                                             material, material, 1, 1, sqp::real(0.01), true);
    sqp::PreparedContactConstraint prepared;
    REQUIRE(sqp::prepare_contact_runtime(bodies, constraint, prepared));
    sqp::solve_contact_constraint(prepared, sqp::real(0.01), true);
    CHECK(velocity(bodies, body).y > 0);
    CHECK(std::isfinite(velocity(bodies, body).y));
}

TEST_CASE("Contact warm start and friction update SoA velocity columns",
          "[contact][warm][friction]") {
    sqp::BodyPool bodies;
    sqp::ShapePool shapes;
    const auto ground =
        add_rigid(bodies, shapes, {0, 0}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto body = add_rigid(bodies, shapes, {0, 1}, {2, -1});
    sqp::PhysicsMaterial material;
    material.base_restitution = sqp::real(0);
    material.static_friction = sqp::real(0.8);
    material.dynamic_friction = sqp::real(0.5);
    auto constraint = prepare_from_materials(bodies, ground, body,
                                             contact({0, sqp::real(0.5)}, {0, 1}, sqp::real(0.01)),
                                             material, material, 1, 1, sqp::real(0.01), false);
    sqp::PreparedContactConstraint prepared;
    REQUIRE(sqp::prepare_contact_runtime(bodies, constraint, prepared));
    sqp::solve_contact_constraint(prepared, 0, false);
    CHECK(std::abs(velocity(bodies, body).x) < sqp::real(2));
    const sqp::real before = velocity(bodies, body).y;
    constraint.points[0].normal_impulse = sqp::real(0.25);
    sqp::warm_start_contact(prepared);
    CHECK(velocity(bodies, body).y == Catch::Approx(before + sqp::real(0.25)));
}

TEST_CASE("Contact friction combines materials by priority", "[contact][friction][material]") {
    sqp::BodyPool bodies;
    sqp::ShapePool shapes;
    const auto ground =
        add_rigid(bodies, shapes, {}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto body = add_rigid(bodies, shapes, {0, 1}, {});
    const auto manifold = contact({0, sqp::real(0.5)}, {0, 1}, sqp::real(0.01));
    sqp::PhysicsMaterial a, b;
    a.static_friction = sqp::real(0.9);
    a.dynamic_friction = sqp::real(0.6);
    b.static_friction = sqp::real(0.4);
    b.dynamic_friction = sqp::real(0.2);

    const auto prepare = [&](const sqp::PhysicsMaterial& first,
                             const sqp::PhysicsMaterial& second) {
        return prepare_from_materials(bodies, ground, body, manifold, first, second, 1, 1,
                                      sqp::real(0.01), false);
    };
    const auto tied = prepare(a, b);
    REQUIRE(tied.count == 1);
    CHECK(tied.static_friction == Catch::Approx(0.4));
    CHECK(tied.dynamic_friction == Catch::Approx(0.2));

    a.friction_priority = 2;
    const auto first_wins = prepare(a, b);
    const auto reversed = prepare(b, a);
    CHECK(first_wins.static_friction == Catch::Approx(0.9));
    CHECK(first_wins.dynamic_friction == Catch::Approx(0.6));
    CHECK(reversed.static_friction == Catch::Approx(first_wins.static_friction));
    CHECK(reversed.dynamic_friction == Catch::Approx(first_wins.dynamic_friction));

    b.friction_priority = 3;
    b.static_friction = sqp::real(0.1);
    b.dynamic_friction = sqp::real(0.3);
    const auto second_wins = prepare(a, b);
    CHECK(second_wins.dynamic_friction == Catch::Approx(0.3));
    CHECK(second_wins.static_friction == Catch::Approx(0.3));
}

TEST_CASE("Contact static and dynamic friction apply bounded tangent impulses",
          "[contact][friction]") {
    struct Case {
        sqp::real initial_x, expected_x, expected_tangent_impulse;
    };
    for (const Case example : {Case{sqp::real(0.25), 0, sqp::real(0.25)},
                               Case{sqp::real(2), sqp::real(1.5), sqp::real(0.5)},
                               Case{sqp::real(-2), sqp::real(-1.5), sqp::real(-0.5)}}) {
        sqp::BodyPool bodies;
        sqp::ShapePool shapes;
        const auto ground =
            add_rigid(bodies, shapes, {}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
        const auto body =
            add_rigid(bodies, shapes, {0, 1}, {example.initial_x, -1},
                      sqp::BodyConfigs::RegisterType::Dynamic,
                      static_cast<sqp::u32>(sqp::BodyConfigs::OverrideFlags::LockRot));
        sqp::PhysicsMaterial material;
        material.base_restitution = 0;
        material.static_friction = sqp::real(0.8);
        material.dynamic_friction = sqp::real(0.5);
        auto constraint = prepare_from_materials(
            bodies, ground, body, contact({0, sqp::real(0.5)}, {0, 1}, sqp::real(0.01)), material,
            material, 1, 1, sqp::real(0.01), false);
        REQUIRE(constraint.count == 1);
        sqp::PreparedContactConstraint prepared;
        REQUIRE(sqp::prepare_contact_runtime(bodies, constraint, prepared));
        sqp::solve_contact_constraint(prepared, 0, false);

        CHECK(constraint.points[0].normal_impulse == Catch::Approx(1));
        CHECK(constraint.points[0].tangent_impulse ==
              Catch::Approx(example.expected_tangent_impulse));
        CHECK(velocity(bodies, body).x == Catch::Approx(example.expected_x));
        CHECK(velocity(bodies, body).y == Catch::Approx(0));
        const auto report = sqp::report_contact_constraint(constraint);
        CHECK(report.points[0].tangent_impulse == Catch::Approx(example.expected_tangent_impulse));
        CHECK(report.tangent_impulse == Catch::Approx(std::abs(example.expected_tangent_impulse)));
        CHECK(report.total_impulse ==
              Catch::Approx(1 + std::abs(example.expected_tangent_impulse)));
    }
}

TEST_CASE("Contact push solves only normals and relaxation solves friction",
          "[contact][friction][bias]") {
    for (bool use_bias : {false, true}) {
        sqp::BodyPool bodies;
        sqp::ShapePool shapes;
        const auto ground =
            add_rigid(bodies, shapes, {}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
        const auto body =
            add_rigid(bodies, shapes, {0, 1}, {2, -1},
                      sqp::BodyConfigs::RegisterType::Dynamic,
                      static_cast<sqp::u32>(sqp::BodyConfigs::OverrideFlags::LockRot));
        sqp::PhysicsMaterial material;
        material.base_restitution = 0;
        material.static_friction = 1;
        material.dynamic_friction = 1;
        auto constraint = prepare_from_materials(
            bodies, ground, body, contact({0, sqp::real(0.5)}, {0, 1}, sqp::real(0.1)), material,
            material, 1, 1, sqp::real(0.01), false);
        sqp::PreparedContactConstraint prepared;
        REQUIRE(sqp::prepare_contact_runtime(bodies, constraint, prepared));
        sqp::solve_contact_constraint(prepared, 0, use_bias);

        const auto normal_impulse = constraint.points[0].normal_impulse;
        CHECK(normal_impulse > 0);
        const auto expected_tangent = use_bias ? sqp::real(0) : normal_impulse;
        CHECK(constraint.points[0].tangent_impulse == Catch::Approx(expected_tangent));
        CHECK(velocity(bodies, body).x == Catch::Approx(2 - expected_tangent));
        CHECK(sqp::report_contact_constraint(constraint).tangent_impulse ==
              Catch::Approx(expected_tangent));
    }
}

TEST_CASE("Contact friction needs a normal impulse in both passes",
          "[contact][friction][bias]") {
    for (bool use_bias : {false, true}) {
        sqp::BodyPool bodies;
        sqp::ShapePool shapes;
        const auto ground =
            add_rigid(bodies, shapes, {}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
        const auto body = add_rigid(
            bodies, shapes, {0, 1}, {2, 0}, sqp::BodyConfigs::RegisterType::Dynamic,
            static_cast<sqp::u32>(sqp::BodyConfigs::OverrideFlags::LockRot));
        sqp::PhysicsMaterial material;
        material.base_restitution = 0;
        material.static_friction = 1;
        material.dynamic_friction = 1;
        auto constraint = prepare_from_materials(
            bodies, ground, body, contact({0, sqp::real(0.5)}, {0, 1}, 0), material,
            material, 1, 1, sqp::real(0.01), false);
        sqp::PreparedContactConstraint prepared;
        REQUIRE(sqp::prepare_contact_runtime(bodies, constraint, prepared));
        sqp::solve_contact_constraint(prepared, 0, use_bias);

        CHECK(constraint.points[0].normal_impulse == 0);
        CHECK(constraint.points[0].tangent_impulse == 0);
        CHECK(velocity(bodies, body).x == Catch::Approx(2));
        CHECK(sqp::report_contact_constraint(constraint).tangent_impulse == 0);
    }
}

TEST_CASE("Contact warm start applies cached tangent impulse", "[contact][friction][warm]") {
    sqp::BodyPool bodies;
    sqp::ShapePool shapes;
    const auto ground =
        add_rigid(bodies, shapes, {}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto body = add_rigid(bodies, shapes, {0, 1}, {}, sqp::BodyConfigs::RegisterType::Dynamic,
                                static_cast<sqp::u32>(sqp::BodyConfigs::OverrideFlags::LockRot));
    sqp::PhysicsMaterial material;
    auto constraint = prepare_from_materials(bodies, ground, body,
                                             contact({0, sqp::real(0.5)}, {0, 1}, sqp::real(0.01)),
                                             material, material, 1, 1, sqp::real(0.01), true);
    REQUIRE(constraint.count == 1);
    constraint.points[0].normal_impulse = sqp::real(0.5);
    constraint.points[0].tangent_impulse = sqp::real(0.2);
    sqp::PreparedContactConstraint prepared;
    REQUIRE(sqp::prepare_contact_runtime(bodies, constraint, prepared));
    sqp::warm_start_contact(prepared);

    CHECK(velocity(bodies, body).x == Catch::Approx(-0.2));
    CHECK(velocity(bodies, body).y == Catch::Approx(0.5));
    CHECK(velocity(bodies, ground) == sqp::Vec2{});
    const auto report = sqp::report_contact_constraint(constraint);
    CHECK(report.normal_impulse == Catch::Approx(0.5));
    CHECK(report.tangent_impulse == Catch::Approx(0.2));
}

TEST_CASE("Two-point contact uses total static support before selecting dynamic friction",
          "[contact][friction][two-point]") {
    sqp::BodyPool bodies;
    sqp::ShapePool shapes;
    const auto ground = add_rigid(bodies, shapes, {}, {},
                                  sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto body = add_rigid(bodies, shapes, {0, 1}, {sqp::real(0.3), 0},
                                sqp::BodyConfigs::RegisterType::Dynamic,
                                static_cast<sqp::u32>(sqp::BodyConfigs::OverrideFlags::LockRot));
    sqp::ContactManifold manifold;
    manifold.count = 2;
    manifold.normal = {0, 1};
    manifold.points[0] = {{sqp::real(-0.5), sqp::real(0.5)}, 0, {}};
    manifold.points[1] = {{sqp::real(0.5), sqp::real(0.5)}, 0, {}};
    sqp::PhysicsMaterial material;
    material.static_friction = 1;
    material.dynamic_friction = sqp::real(0.1);
    auto constraint = prepare_from_materials(bodies, ground, body, manifold, material,
                                              material, 1, 1, sqp::real(0.1), true);
    constraint.points[0].normal_impulse = sqp::real(0.1);
    constraint.points[1].normal_impulse = sqp::real(0.9);
    constraint.points[0].tangent_impulse = sqp::real(0.1);
    sqp::PreparedContactConstraint prepared;
    REQUIRE(sqp::prepare_contact_runtime(bodies, constraint, prepared));
    sqp::solve_contact_constraint(prepared, 0, false);
    CHECK(constraint.points[0].tangent_impulse == Catch::Approx(sqp::real(0.1)));
    CHECK(constraint.points[1].tangent_impulse > sqp::real(0.01));
    CHECK(velocity(bodies, body).x == Catch::Approx(0).margin(1e-5));
}

TEST_CASE("True contact sliding uses dynamic friction", "[contact][friction]") {
    sqp::BodyPool bodies;
    sqp::ShapePool shapes;
    const auto ground = add_rigid(bodies, shapes, {}, {},
                                  sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto body = add_rigid(bodies, shapes, {0, 1}, {3, 0},
                                sqp::BodyConfigs::RegisterType::Dynamic,
                                static_cast<sqp::u32>(sqp::BodyConfigs::OverrideFlags::LockRot));
    sqp::PhysicsMaterial material;
    material.static_friction = 1;
    material.dynamic_friction = sqp::real(0.5);
    auto constraint = prepare_from_materials(bodies, ground, body,
                                              contact({0, sqp::real(0.5)}, {0, 1}, 0),
                                              material, material, 1, 1, sqp::real(0.1), true);
    auto& point = constraint.points[0];
    point.normal_impulse = 1;
    sqp::BodyState state;
    REQUIRE(bodies.get_state(body, state));
    state.position.x = sqp::real(0.01);
    REQUIRE(bodies.set_state(body, state));
    sqp::PreparedContactConstraint prepared;
    REQUIRE(sqp::prepare_contact_runtime(bodies, constraint, prepared));
    sqp::refresh_contact_constraint_geometry(prepared);
    sqp::solve_contact_constraint(prepared, 0, false);
    CHECK(point.tangent_impulse == Catch::Approx(sqp::real(0.5)));
    CHECK(velocity(bodies, body).x == Catch::Approx(sqp::real(2.5)));
}

TEST_CASE("Static friction has no positional anchor correction", "[contact][friction]") {
    for (const sqp::real shift : {sqp::real(0.01), sqp::real(0.5)}) {
        sqp::BodyPool bodies;
        sqp::ShapePool shapes;
        const auto ground =
            add_rigid(bodies, shapes, {}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
        const auto body = add_rigid(
            bodies, shapes, {0, 1}, {}, sqp::BodyConfigs::RegisterType::Dynamic,
            static_cast<sqp::u32>(sqp::BodyConfigs::OverrideFlags::LockRot));
        sqp::PhysicsMaterial material;
        material.static_friction = 1;
        material.dynamic_friction = sqp::real(0.5);
        auto constraint = prepare_from_materials(
            bodies, ground, body, contact({0, sqp::real(0.5)}, {0, 1}, 0), material,
            material, 1, 1, sqp::real(0.1), true);
        auto& point = constraint.points[0];
        point.normal_impulse = 1;
        sqp::BodyState state;
        REQUIRE(bodies.get_state(body, state));
        state.position.x = shift;
        REQUIRE(bodies.set_state(body, state));
        sqp::PreparedContactConstraint prepared;
        REQUIRE(sqp::prepare_contact_runtime(bodies, constraint, prepared));
        sqp::refresh_contact_constraint_geometry(prepared);
        sqp::solve_contact_constraint(prepared, 0, false);

        CHECK(velocity(bodies, body).x == Catch::Approx(0));
        CHECK(point.tangent_impulse == Catch::Approx(0));
    }
}

TEST_CASE("Contact friction requires normal support", "[contact][friction]") {
    sqp::BodyPool bodies;
    sqp::ShapePool shapes;
    const auto ground =
        add_rigid(bodies, shapes, {}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto body = add_rigid(bodies, shapes, {0, 1}, {});
    sqp::PhysicsMaterial material;
    material.static_friction = 1;
    auto constraint = prepare_from_materials(
        bodies, ground, body, contact({0, sqp::real(0.5)}, {0, 1}, 0), material,
        material, 1, 1, sqp::real(0.1), true);
    sqp::BodyState state;
    REQUIRE(bodies.get_state(body, state));
    state.position.x = sqp::real(0.01);
    REQUIRE(bodies.set_state(body, state));
    sqp::PreparedContactConstraint prepared;
    REQUIRE(sqp::prepare_contact_runtime(bodies, constraint, prepared));
    sqp::refresh_contact_constraint_geometry(prepared);
    sqp::solve_contact_constraint(prepared, 0, false);
    CHECK(velocity(bodies, body).x == Catch::Approx(0));
    CHECK(constraint.points[0].tangent_impulse == 0);
}

TEST_CASE("Contact supports lite bodies and honors locked scalar axes", "[contact][lite][lock]") {
    sqp::BodyPool bodies;
    sqp::ShapePool shapes;
    const auto locked =
        add_lite(bodies, shapes, {0, 0}, {},
                 static_cast<sqp::u32>(sqp::BodyConfigs::OverrideFlags::LockLinearX));
    const auto moving = add_lite(bodies, shapes, {2, 0}, {-1, 0});
    sqp::PhysicsMaterial material;
    auto constraint =
        prepare_from_materials(bodies, locked, moving, contact({1, 0}, {1, 0}, sqp::real(0.1)),
                               material, material, 1, 1, sqp::real(0.01), false);
    REQUIRE(constraint.count == 1);
    CHECK(constraint.inverse_mass_a.x == 0);
    CHECK(constraint.inverse_inertia_a == 0);
    sqp::PreparedContactConstraint prepared;
    REQUIRE(sqp::prepare_contact_runtime(bodies, constraint, prepared));
    sqp::solve_contact_constraint(prepared, 0, false);
    CHECK(velocity(bodies, locked).x == 0);
    CHECK(std::isfinite(velocity(bodies, moving).x));
}

TEST_CASE("Contact handles rigid-lite two-point contacts and zero effective mass",
          "[contact][lite][two-point]") {
    sqp::BodyPool bodies;
    sqp::ShapePool shapes;
    const auto rigid =
        add_rigid(bodies, shapes, {0, 0}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto lite = add_lite(bodies, shapes, {0, 1}, {0, -1});
    sqp::ContactManifold manifold;
    manifold.count = 2;
    manifold.normal = {0, 1};
    manifold.points[0] = {{sqp::real(-0.5), sqp::real(0.5)}, sqp::real(0.05), {}};
    manifold.points[1] = {{sqp::real(0.5), sqp::real(0.5)}, sqp::real(0.05), {}};
    sqp::PhysicsMaterial material;
    material.base_restitution = sqp::real(0);
    auto constraint = prepare_from_materials(bodies, rigid, lite, manifold, material, material, 1,
                                             1, sqp::real(0.01), false);
    REQUIRE(constraint.count == 2);
    sqp::PreparedContactConstraint prepared;
    REQUIRE(sqp::prepare_contact_runtime(bodies, constraint, prepared));
    sqp::solve_contact_constraint(prepared, 0, false);
    CHECK(std::isfinite(velocity(bodies, lite).y));
    CHECK(velocity(bodies, lite).y >= sqp::real(-1e-5));

    const auto other_static =
        add_rigid(bodies, shapes, {0, 2}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    auto immovable =
        prepare_from_materials(bodies, rigid, other_static, contact({0, 1}, {0, 1}, sqp::real(0.1)),
                               material, material, 1, 1, sqp::real(0.01), false);
    REQUIRE(immovable.count == 1);
    CHECK(immovable.points[0].normal_mass == 0);
    sqp::PreparedContactConstraint prepared_immovable;
    REQUIRE(sqp::prepare_contact_runtime(bodies, immovable, prepared_immovable));
    sqp::solve_contact_constraint(prepared_immovable, 0, false);
    CHECK(velocity(bodies, rigid) == sqp::Vec2{});
    CHECK(velocity(bodies, other_static) == sqp::Vec2{});
}

TEST_CASE("Contact rejects invalid inputs without changing bodies", "[contact][validation]") {
    sqp::BodyPool bodies;
    sqp::ShapePool shapes;
    const auto body = add_rigid(bodies, shapes, {}, {3, 4});
    sqp::PhysicsMaterial material;
    sqp::ContactManifold empty;
    auto invalid = prepare_from_materials(bodies, body, {}, empty, material, material, 1, 1,
                                          sqp::real(0.01), false);
    CHECK(invalid.count == 0);
    auto zero_step = prepare_from_materials(bodies, body, body, contact({}, {0, 1}, 1), material,
                                            material, 1, 1, 0, false);
    CHECK(zero_step.count == 0);
    sqp::PreparedContactConstraint prepared;
    CHECK_FALSE(sqp::prepare_contact_runtime(bodies, invalid, prepared));
    sqp::solve_contact_constraint(prepared, 0, false);
    CHECK(velocity(bodies, body) == sqp::Vec2(3, 4));

    sqp::BodyState sleeping;
    REQUIRE(bodies.get_state(body, sleeping));
    REQUIRE(sqp::prepare_contact_runtime(bodies, zero_step, prepared));
    sqp::BodyState after;
    REQUIRE(bodies.get_state(body, after));
    CHECK(after.flags == sleeping.flags);
}

TEST_CASE("Contact dynamic restitution uses the lesser shape impact scale",
          "[contact][restitution]") {
    sqp::BodyPool bodies;
    sqp::ShapePool shapes;
    const auto ground =
        add_rigid(bodies, shapes, {0, 0}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto body = add_rigid(bodies, shapes, {0, 1}, {0, -2});
    sqp::PhysicsMaterial material;
    material.base_restitution = sqp::real(1);
    material.restitution_impact_threshold = sqp::real(1);
    material.restitution_impact_exponent = sqp::real(1);

    const auto manifold = contact({0, sqp::real(0.5)}, {0, 1}, sqp::real(0.1));
    const auto first = prepare_from_materials(bodies, ground, body, manifold, material, material, 2,
                                              4, sqp::real(0.01), false);
    const auto reversed = prepare_from_materials(bodies, ground, body, manifold, material, material,
                                                 4, 2, sqp::real(0.01), false);
    const auto both_high = prepare_from_materials(bodies, ground, body, manifold, material,
                                                  material, 4, 4, sqp::real(0.01), false);
    REQUIRE(first.count == 1);
    CHECK(first.restitution == Catch::Approx(reversed.restitution));
    CHECK(first.restitution > both_high.restitution);
}

TEST_CASE("Separated speculative contacts do not apply restitution", "[contact][restitution]") {
    sqp::BodyPool bodies;
    sqp::ShapePool shapes;
    const auto ground =
        add_rigid(bodies, shapes, {0, 0}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto body = add_rigid(bodies, shapes, {0, 1}, {0, -10});
    sqp::PhysicsMaterial material;
    material.base_restitution = sqp::real(1);
    const auto separated = prepare_from_materials(
        bodies, ground, body, contact({0, sqp::real(0.5)}, {0, 1}, sqp::real(-0.1)), material,
        material, 1, 1, sqp::real(0.01), false);
    const auto touching =
        prepare_from_materials(bodies, ground, body, contact({0, sqp::real(0.5)}, {0, 1}, 0),
                               material, material, 1, 1, sqp::real(0.01), false);
    REQUIRE(separated.count == 1);
    REQUIRE(touching.count == 1);
    CHECK(separated.points[0].restitution_velocity == 0);
    CHECK(touching.points[0].restitution_velocity > 0);
}
