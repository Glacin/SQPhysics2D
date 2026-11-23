#include <cmath>
#include <limits>

#include "ccd/convex_toi.h"
#include "ccd/shape_ccd.h"
#include "narrow_phase/sat.h"
#include "sqphysics2d/detail/data_storage.h"
#include <catch_amalgamated.hpp>

namespace sqp {

// Keeps each numerical scenario in a single pool while allowing the existing
// test cases to name shapes independently.
struct TestShape {
    static ShapePool& pool() {
        static ShapePool instance;
        return instance;
    }

    ShapeHandle handle;
    Vec2 local_center{std::numeric_limits<real>::quiet_NaN(),
                      std::numeric_limits<real>::quiet_NaN()};
    real inner_bound_radius = -1;
    real bound_radius_about_centroid = -1;

    void set_geometry(const Circle& circle, Vec2 offset = {}) {
        if (is_valid(handle))
            pool().remove(handle);
        CircleShapeDef def;
        def.radius = circle.radius;
        def.local_offset = circle.center + offset;
        handle = pool().add_circle(def);
    }
    void set_geometry(const ConvexPoly& polygon, Vec2 offset = {}) {
        if (is_valid(handle))
            pool().remove(handle);
        PolygonShapeDef def;
        def.local_offset = offset + polygon.center_offset();
        for (Vec2 vertex : polygon.local)
            def.vertices.push_back(vertex);
        handle = pool().add_polygon(def);
    }
    operator detail::ShapeQueryView() const {
        detail::ShapeQueryView view(pool(), handle);
        if (math::isfinite(local_center.x))
            view.local_center = local_center;
        if (inner_bound_radius >= 0)
            view.inner_bound_radius = inner_bound_radius;
        if (bound_radius_about_centroid >= 0)
            view.bound_radius_about_centroid = bound_radius_about_centroid;
        return view;
    }
};

} // namespace sqp

namespace {

bool finite_vec(const sqp::Vec2& value) {
    return sqp::math::isfinite(value.x) && sqp::math::isfinite(value.y);
}

sqp::TestShape unit_box_shape() {
    sqp::TestShape shape;
    shape.set_geometry(sqp::ConvexPoly({{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}));
    shape.inner_bound_radius = 1;
    shape.bound_radius_about_centroid = sqp::real(1.5);
    return shape;
}

} // namespace

TEST_CASE("contact kinematics combines scalar linear and angular velocities",
          "[collision][ccd][kinematics]") {
    sqp::TestShape a, b;
    a.set_geometry(sqp::Circle({1, 0}, 1));
    b.set_geometry(sqp::Circle({1, 0}, 1));
    const sqp::Vec2 position_a(10, 20), position_b(12, 20);
    const sqp::real rotation = sqp::math::PI * sqp::real(0.5);
    sqp::Vec2 velocity_a, velocity_b;
    sqp::real angular_a = 0, angular_b = 0;
    sqp::real expected_speed = 0;
    SECTION("linear only") {
        velocity_a = {2, 0};
        velocity_b = {-1, 0};
        expected_speed = -3;
    }
    SECTION("angular only") {
        angular_a = 3;
        angular_b = -2;
        expected_speed = 5;
    }
    SECTION("combined") {
        velocity_a = {2, 0};
        velocity_b = {-1, 0};
        angular_a = 3;
        angular_b = -2;
        expected_speed = 2;
    }
    const auto result =
        sqp::ccd_contact_kinematics(a, position_a, rotation, velocity_a, angular_a, b, position_b,
                                    rotation, velocity_b, angular_b, {-1, 0});
    CHECK(result.point.x == Catch::Approx(11).margin(1e-5));
    CHECK(result.point.y == Catch::Approx(21).margin(1e-5));
    CHECK(result.normal_ab.x == Catch::Approx(1));
    CHECK(result.normal_ab.y == Catch::Approx(0));
    CHECK(result.normal_speed == Catch::Approx(expected_speed).margin(1e-5));
}

TEST_CASE("A zero-time polygon sweep reports the touching wall normal",
          "[collision][ccd][shapecast][regression]") {
    sqp::TestShape wall, body;
    wall.set_geometry(sqp::ConvexPoly({{-129, -0.5f}, {129, -0.5f}, {129, 0.5f}, {-129, 0.5f}}));
    body.set_geometry(
        sqp::ConvexPoly({{-0.25f, -0.25f}, {0.25f, -0.25f}, {0.25f, 0.25f}, {-0.25f, 0.25f}}));
    const sqp::Vec2 start{-114, -127.75f}, end{-115, -129}, wall_position{0, -128.5f};
    bool sweeps_valid = true;
    for (const auto rotation : {sqp::real(0), sqp::real(0.1)}) {
        sqp::real fraction = 1;
        sqp::Vec2 normal;
        sweeps_valid &=
            sqp::sweep_shape_pair_toi(body, start, 0, end, rotation, wall, wall_position, 0,
                                      wall_position, 0, true, 0, fraction, normal) &&
            fraction == 0 && normal.x == Catch::Approx(0).margin(1e-6) &&
            normal.y == Catch::Approx(1).margin(1e-6) && sqp::math::dot(end - start, normal) < 0 &&
            sqp::math::dot(sqp::Vec2{1, 0}, normal) == 0 &&
            sqp::math::dot(sqp::Vec2{0, 1}, normal) > 0;
    }
    CHECK(sweeps_valid);
}

TEST_CASE("scalar CCD entries preserve statistics", "[collision][ccd][statistics]") {
    sqp::TestShape a, b;
    a.set_geometry(sqp::Circle({3, 0}, sqp::real(0.5)));
    b.set_geometry(sqp::Circle({0, 3}, sqp::real(0.5)));
    const sqp::Vec2 position(10, 20);
    const sqp::real rotation1 = sqp::math::PI * sqp::real(0.5);
    const sqp::real penetration = GENERATE(sqp::real(0), sqp::real(1));
    const auto reference = sqp::convex_time_of_impact(a, position, 0, position, rotation1, b,
                                                      position, 0, position, 0, penetration);
    const bool conservative_failure = reference.state == sqp::ConvexToiState::Failed;
    REQUIRE((reference.state == sqp::ConvexToiState::Hit || conservative_failure));
    REQUIRE(reference.fraction < 1);
    const bool enabled = GENERATE(false, true);
    const bool cached = GENERATE(false, true);
    sqp::real toi = -1;
    sqp::Vec2 normal;
    sqp::u32 distance_iterations = 999, root_iterations = 999;
    bool failed = true;
    bool hit;
    if (cached) {
        sqp::PreparedPolygon prepared_a, prepared_b;
        hit = sqp::sweep_shape_pair_toi(
            a, position, 0, position, rotation1, b, position, 0, position, 0, enabled, penetration,
            toi, normal, &distance_iterations, &root_iterations, &failed, prepared_a, prepared_b);
    } else {
        hit = sqp::sweep_shape_pair_toi(a, position, 0, position, rotation1, b, position, 0,
                                        position, 0, enabled, penetration, toi, normal,
                                        &distance_iterations, &root_iterations, &failed);
    }
    CHECK(hit == enabled);
    CHECK(failed == (enabled && conservative_failure));
    CHECK(toi == (enabled ? reference.fraction : sqp::real(1)));
    CHECK(normal == (enabled ? reference.normal : sqp::Vec2(1, 0)));
    CHECK(distance_iterations == (enabled ? reference.distance_iterations : 0));
    CHECK(root_iterations == (enabled ? reference.root_iterations : 0));
}

TEST_CASE("offset geometry preserves CCD for all shape pairs", "[collision][ccd][centroid]") {
    const bool circle_a = GENERATE(false, true);
    const bool circle_b = GENERATE(false, true);
    CAPTURE(circle_a, circle_b);
    const auto make_shape = [](bool circle) {
        sqp::TestShape shape;
        if (circle)
            shape.set_geometry(sqp::Circle({2, 1}, 1), {1, -1});
        else
            shape.set_geometry(sqp::ConvexPoly({{1, 0}, {3, 0}, {3, 2}, {1, 2}}), {1, -1});
        return shape;
    };
    const auto a = make_shape(circle_a);
    const auto b = make_shape(circle_b);

    sqp::Vec2 position_a0;
    sqp::real rotation_a0 = 0;
    sqp::Vec2 position_a1;
    sqp::real rotation_a1 = 0;

    sqp::Vec2 position_b0;
    sqp::real rotation_b0 = 0;
    sqp::Vec2 position_b1;
    sqp::real rotation_b1 = 0;

    position_a0 = {10, 20};
    position_a1 = {20, 20};
    position_b0 = position_b1 = {16, 26};
    rotation_a0 = rotation_a1 = sqp::math::PI * sqp::real(0.5);
    rotation_b0 = rotation_b1 = -sqp::math::PI * sqp::real(0.5);

    bool enabled = true;
    sqp::real target_penetration_m = 1;
    target_penetration_m = 0;
    sqp::real toi = 1;
    sqp::Vec2 normal;
    REQUIRE(sqp::sweep_shape_pair_toi(a, position_a0, rotation_a0, position_a1, rotation_a1, b,
                                      position_b0, rotation_b0, position_b1, rotation_b1, enabled,
                                      target_penetration_m, toi, normal));
    CHECK(toi == Catch::Approx(0.4).margin(1e-4));
    CHECK(normal.x == Catch::Approx(-1).margin(1e-4));
    // Confirm the shape cast agrees with the distance-proxy TOI result.
    const auto result = sqp::convex_time_of_impact(a, position_a0, rotation_a0, position_a1,
                                                   rotation_a1, b, position_b0, rotation_b0,
                                                   position_b1, rotation_b1, target_penetration_m);
    REQUIRE(result.state == sqp::ConvexToiState::Hit);

    auto centered_a = a;
    auto centered_b = b;

    centered_a.local_center = centered_b.local_center = {};
    auto centered_pose_a_position0 = position_a0;
    auto centered_pose_a_rotation0 = rotation_a0;

    auto centered_pose_a_position1 = position_a1;
    auto centered_pose_a_rotation1 = rotation_a1;

    auto centered_pose_b_position0 = position_b0;
    auto centered_pose_b_rotation0 = rotation_b0;

    auto centered_pose_b_position1 = position_b1;
    auto centered_pose_b_rotation1 = rotation_b1;

    centered_pose_a_position0 += sqp::Vec2(0, 3);
    centered_pose_a_position1 += sqp::Vec2(0, 3);
    centered_pose_b_position0 -= sqp::Vec2(0, 3);
    centered_pose_b_position1 -= sqp::Vec2(0, 3);
    const auto reference = sqp::convex_time_of_impact(
        centered_a, centered_pose_a_position0, centered_pose_a_rotation0, centered_pose_a_position1,
        centered_pose_a_rotation1, centered_b, centered_pose_b_position0, centered_pose_b_rotation0,
        centered_pose_b_position1, centered_pose_b_rotation1, target_penetration_m);
    CHECK(result.state == reference.state);
    CHECK(result.fraction == Catch::Approx(reference.fraction).margin(1e-4));
}

TEST_CASE("offset circle rotation detects impact and reentry", "[collision][ccd][centroid]") {
    sqp::TestShape a, b;
    a.set_geometry(sqp::Circle({3, 0}, sqp::real(0.5)));
    b.set_geometry(sqp::Circle({0, 3}, sqp::real(0.5)));
    sqp::Vec2 position_a0;
    sqp::real rotation_a0 = 0;
    sqp::Vec2 position_a1;
    sqp::real rotation_a1 = 0;
    sqp::Vec2 position_b0;
    sqp::real rotation_b0 = 0;
    sqp::Vec2 position_b1;
    sqp::real rotation_b1 = 0;
    position_a0 = position_a1 = position_b0 = position_b1 = {10, 20};
    rotation_a1 = sqp::math::PI * sqp::real(0.5);
    bool enabled = true;
    sqp::real target_penetration_m = 0;
    sqp::real toi = 1;
    sqp::Vec2 normal;
    REQUIRE(sqp::sweep_shape_pair_toi(a, position_a0, rotation_a0, position_a1, rotation_a1, b,
                                      position_b0, rotation_b0, position_b1, rotation_b1, enabled,
                                      target_penetration_m, toi, normal));
    CHECK(toi > sqp::real(0.7));
    CHECK(toi < sqp::real(0.9));
    rotation_a0 = sqp::math::PI * sqp::real(0.5);
    rotation_a1 = sqp::math::PI * sqp::real(2.5);
    REQUIRE(sqp::sweep_shape_pair_reentry_toi(a, position_a0, rotation_a0, position_a1, rotation_a1,
                                              b, position_b0, rotation_b0, position_b1, rotation_b1,
                                              enabled, target_penetration_m, toi, normal));
    CHECK(toi > sqp::real(0.9));
    CHECK(toi < sqp::real(1));
}

TEST_CASE("rotating circle reaches polygon wall after its initial gap", "[collision][ccd][toi]") {
    sqp::TestShape circle, wall;
    circle.set_geometry(sqp::Circle({1, 0}, sqp::real(0.5)));
    wall.set_geometry(sqp::ConvexPoly({{-0.1, -2}, {0.1, -2}, {0.1, 2}, {-0.1, 2}}));
    const sqp::Vec2 position(0, 0);
    const sqp::real quarter_turn = sqp::math::PI * sqp::real(0.5);
    const auto contact = sqp::convex_time_of_impact(circle, position, 0, position, quarter_turn,
                                                    wall, position, 0, position, 0, 0);
    const auto allowed = sqp::convex_time_of_impact(circle, position, 0, position, quarter_turn,
                                                    wall, position, 0, position, 0, sqp::real(0.1));
    REQUIRE(contact.state == sqp::ConvexToiState::Hit);
    REQUIRE(allowed.state == sqp::ConvexToiState::Hit);
    CHECK(contact.fraction > sqp::real(0));
    CHECK(contact.fraction == Catch::Approx(std::acos(0.6) / quarter_turn).margin(0.002));
    CHECK(allowed.fraction == Catch::Approx(std::acos(0.5) / quarter_turn).margin(0.002));
    CHECK(allowed.fraction > contact.fraction);
}

TEST_CASE("scalar TOI keeps the first hit of a full rotation", "[collision][ccd][toi]") {
    sqp::TestShape circle, wall;
    circle.set_geometry(sqp::Circle({1, 0}, sqp::real(0.2)));
    wall.set_geometry(sqp::ConvexPoly({{-0.1, -2}, {0.1, -2}, {0.1, 2}, {-0.1, 2}}));
    const sqp::Vec2 position(0, 0);
    const sqp::real full_turn = sqp::math::PI * sqp::real(2);
    const sqp::real first_contact = std::acos(sqp::real(0.3)) / full_turn;
    const auto first = sqp::convex_time_of_impact(circle, position, 0, position, full_turn, wall,
                                                  position, 0, position, 0, 0);
    REQUIRE(first.state == sqp::ConvexToiState::Hit);
    CHECK(first.fraction == Catch::Approx(first_contact).margin(0.002));
    CHECK(first.normal.x > sqp::real(0.9));
    CHECK(first.distance_iterations > 0);

    // The sweep is separated again at its end and has a second entry after
    // crossing the wall. A bounded query must find that later entry.
    const auto reentry =
        sqp::convex_time_of_impact(circle, position, 0, position, full_turn, wall, position, 0,
                                   position, 0, sqp::real(0), sqp::real(0.4), sqp::real(1));
    REQUIRE(reentry.state == sqp::ConvexToiState::Hit);
    CHECK(reentry.fraction ==
          Catch::Approx(sqp::real(1) - std::acos(sqp::real(-0.3)) / full_turn).margin(0.002));
    CHECK(reentry.normal.x < sqp::real(-0.9));
}

TEST_CASE("scalar TOI respects overlap miss and bounded sweeps", "[collision][ccd][toi]") {
    sqp::TestShape circle, wall;
    circle.set_geometry(sqp::Circle({1, 0}, sqp::real(0.2)));
    wall.set_geometry(sqp::ConvexPoly({{-0.1, -2}, {0.1, -2}, {0.1, 2}, {-0.1, 2}}));
    const sqp::Vec2 position(0, 0);
    const sqp::real full_turn = sqp::math::PI * sqp::real(2);
    const auto overlap =
        sqp::convex_time_of_impact(circle, position, 0, position, full_turn, wall, position, 0,
                                   position, 0, sqp::real(0), sqp::real(0.25));
    CHECK(overlap.state == sqp::ConvexToiState::InitialOverlap);
    CHECK(overlap.fraction == Catch::Approx(0.25));
    const auto bounded =
        sqp::convex_time_of_impact(circle, position, 0, position, full_turn, wall, position, 0,
                                   position, 0, sqp::real(0), sqp::real(0), sqp::real(0.1));
    CHECK(bounded.state == sqp::ConvexToiState::Separated);
    const auto miss = sqp::convex_time_of_impact(circle, position, 0, position, full_turn, wall,
                                                 {3, 0}, 0, {3, 0}, 0, 0);
    CHECK(miss.state == sqp::ConvexToiState::Separated);
}

TEST_CASE("empty and invalid variants cannot collide", "[collision][ccd][centroid]") {
    sqp::TestShape invalid;
    SECTION("monostate") {}
    SECTION("polygon") {
        invalid.set_geometry(sqp::ConvexPoly());
    }
    SECTION("circle") {
        invalid.set_geometry(sqp::Circle({2, 3}, 0));
    }
    const auto valid = unit_box_shape();
    sqp::BodyState state;
    CHECK_FALSE(sqp::convex_signed_distance(invalid, state.position, state.rotation, valid,
                                            state.position, state.rotation)
                    .valid);
    sqp::Vec2 pose_position0;
    sqp::real pose_rotation0 = 0;
    sqp::Vec2 pose_position1;
    sqp::real pose_rotation1 = 0;
    pose_rotation1 = 1;
    bool enabled = true;
    sqp::real target_penetration_m = 1;
    sqp::real toi = 1;
    sqp::Vec2 normal;
    CHECK_FALSE(sqp::sweep_shape_pair_toi(invalid, pose_position0, pose_rotation0, pose_position1,
                                          pose_rotation1, valid, pose_position0, pose_rotation0,
                                          pose_position1, pose_rotation1, enabled,
                                          target_penetration_m, toi, normal));
    CHECK_FALSE(sqp::sweep_shape_pair_reentry_toi(
        invalid, pose_position0, pose_rotation0, pose_position1, pose_rotation1, valid,
        pose_position0, pose_rotation0, pose_position1, pose_rotation1, enabled,
        target_penetration_m, toi, normal));
}

TEST_CASE("convex distance exposes separated surface witnesses", "[collision][distance]") {
    const sqp::TestShape a = unit_box_shape();
    const sqp::TestShape b = unit_box_shape();
    sqp::BodyState state_a, state_b;
    state_b.position = {3, 0};
    const auto distance = sqp::convex_signed_distance(a, state_a.position, state_a.rotation, b,
                                                      state_b.position, state_b.rotation);
    REQUIRE(distance.valid);
    CHECK_FALSE(distance.overlapped);
    CHECK(distance.separation == Catch::Approx(1));
    CHECK(distance.normal.x == Catch::Approx(1));
    CHECK(distance.normal.y == Catch::Approx(0));
    CHECK(distance.point_a.x == Catch::Approx(1));
    CHECK(distance.point_b.x == Catch::Approx(2));
}

TEST_CASE("polygon translation shape cast returns earliest impact", "[collision][ccd]") {
    const sqp::TestShape a = unit_box_shape();
    const sqp::TestShape b = unit_box_shape();

    sqp::Vec2 position_a0;
    sqp::real rotation_a0 = 0;
    sqp::Vec2 position_a1;
    sqp::real rotation_a1 = 0;
    position_a0 = {0, 0};
    position_a1 = {10, 0};
    sqp::Vec2 position_b0;
    sqp::real rotation_b0 = 0;
    sqp::Vec2 position_b1;
    sqp::real rotation_b1 = 0;
    position_b0 = position_b1 = {6, 0};
    bool enabled = true;
    sqp::real target_penetration_m = 1;
    target_penetration_m = 0;
    sqp::real toi = 1;
    sqp::Vec2 normal;
    REQUIRE(sqp::sweep_shape_pair_toi(a, position_a0, rotation_a0, position_a1, rotation_a1, b,
                                      position_b0, rotation_b0, position_b1, rotation_b1, enabled,
                                      target_penetration_m, toi, normal));
    CHECK(toi == Catch::Approx(0.4f));
    CHECK(normal.x == Catch::Approx(-1));
}

TEST_CASE("shape cast gives the correct normal for offset polygons", "[collision][ccd]") {
    sqp::TestShape moving, target;
    moving.set_geometry(sqp::ConvexPoly({{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}), {-10, 0});
    target.set_geometry(sqp::ConvexPoly({{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}));
    sqp::real toi = 1;
    sqp::Vec2 normal;
    REQUIRE(sqp::sweep_shape_pair_toi(moving, {10, 0}, 0, {20, 0}, 0, target, {6, 0}, 0, {6, 0}, 0,
                                      true, 0, toi, normal));
    CHECK(toi == Catch::Approx(sqp::real(0.4)).margin(1e-4));
    CHECK(normal.x == Catch::Approx(-1).margin(1e-4));
}

TEST_CASE("shapecast handles rotating polygon sweep", "[collision][ccd][shapecast]") {
    sqp::TestShape blade;
    sqp::TestShape target;
    blade.set_geometry(sqp::ConvexPoly({{-2, -0.1f}, {2, -0.1f}, {2, 0.1f}, {-2, 0.1f}}));
    target.set_geometry(
        sqp::ConvexPoly({{-0.5f, -0.5f}, {0.5f, -0.5f}, {0.5f, 0.5f}, {-0.5f, 0.5f}}));
    blade.inner_bound_radius = sqp::real(0.1);
    target.inner_bound_radius = sqp::real(0.5);

    sqp::Vec2 position_a0;
    sqp::real rotation_a0 = 0;
    sqp::Vec2 position_a1;
    sqp::real rotation_a1 = 0;
    position_a0 = {0, 0};
    position_a1 = {0, 0};
    rotation_a0 = 0;
    rotation_a1 = sqp::math::PI * sqp::real(0.5);
    sqp::Vec2 position_b0;
    sqp::real rotation_b0 = 0;
    sqp::Vec2 position_b1;
    sqp::real rotation_b1 = 0;
    position_b0 = position_b1 = {0, sqp::real(1.5)};
    bool enabled = true;
    sqp::real target_penetration_m = 1;
    target_penetration_m = 0;
    sqp::real toi = 1;
    sqp::Vec2 normal;

    REQUIRE(sqp::sweep_shape_pair_toi(blade, position_a0, rotation_a0, position_a1, rotation_a1,
                                      target, position_b0, rotation_b0, position_b1, rotation_b1,
                                      enabled, target_penetration_m, toi, normal));
    CHECK(toi >= sqp::real(0));
    CHECK(toi <= sqp::real(1));
    CHECK(finite_vec(normal));
}

TEST_CASE("shapecast preserves rotations larger than pi", "[collision][ccd][shapecast][angular]") {
    sqp::TestShape blade;
    blade.set_geometry(sqp::ConvexPoly({{-2, -0.1f}, {2, -0.1f}, {2, 0.1f}, {-2, 0.1f}}));
    blade.inner_bound_radius = sqp::real(0.1);
    blade.bound_radius_about_centroid = sqp::real(2.01);

    sqp::TestShape target = unit_box_shape();
    target.set_geometry(
        sqp::ConvexPoly({{-0.25f, -0.25f}, {0.25f, -0.25f}, {0.25f, 0.25f}, {-0.25f, 0.25f}}));
    target.inner_bound_radius = sqp::real(0.25);
    target.bound_radius_about_centroid = sqp::real(0.36);

    sqp::Vec2 position_a0;
    sqp::real rotation_a0 = 0;
    sqp::Vec2 position_a1;
    sqp::real rotation_a1 = 0;
    position_a0 = position_a1 = {0, 0};
    rotation_a0 = 0;
    rotation_a1 = sqp::math::PI * sqp::real(2.5);
    sqp::Vec2 position_b0;
    sqp::real rotation_b0 = 0;
    sqp::Vec2 position_b1;
    sqp::real rotation_b1 = 0;
    position_b0 = position_b1 = {0, sqp::real(1.5)};
    bool enabled = true;
    sqp::real target_penetration_m = 1;
    sqp::real toi = 1;
    sqp::Vec2 normal;

    REQUIRE(sqp::sweep_shape_pair_toi(blade, position_a0, rotation_a0, position_a1, rotation_a1,
                                      target, position_b0, rotation_b0, position_b1, rotation_b1,
                                      enabled, target_penetration_m, toi, normal));
    CHECK(toi > sqp::real(0));
    CHECK(toi < sqp::real(0.25));
    CHECK(finite_vec(normal));
}

TEST_CASE("shapecast finds rotation reentry after initial separation",
          "[collision][ccd][shapecast][angular][reentry]") {
    sqp::TestShape blade;
    blade.set_geometry(sqp::ConvexPoly({{-2, -0.1f}, {2, -0.1f}, {2, 0.1f}, {-2, 0.1f}}));
    blade.inner_bound_radius = sqp::real(0.1);
    blade.bound_radius_about_centroid = sqp::real(2.01);
    sqp::TestShape target = unit_box_shape();

    sqp::Vec2 position_a0;
    sqp::real rotation_a0 = 0;
    sqp::Vec2 position_a1;
    sqp::real rotation_a1 = 0;
    position_a0 = position_a1 = {0, 0};
    rotation_a0 = sqp::math::PI * sqp::real(0.5);
    rotation_a1 = sqp::math::PI * sqp::real(2.5);
    sqp::Vec2 position_b0;
    sqp::real rotation_b0 = 0;
    sqp::Vec2 position_b1;
    sqp::real rotation_b1 = 0;
    position_b0 = position_b1 = {0, sqp::real(1.5)};
    bool enabled = true;
    sqp::real target_penetration_m = 1;
    sqp::real toi = 1;
    sqp::Vec2 normal;

    REQUIRE(sqp::sweep_shape_pair_reentry_toi(
        blade, position_a0, rotation_a0, position_a1, rotation_a1, target, position_b0, rotation_b0,
        position_b1, rotation_b1, enabled, target_penetration_m, toi, normal));
    CHECK(toi > sqp::real(0.25));
    CHECK(toi < sqp::real(0.75));
    CHECK(finite_vec(normal));
}

TEST_CASE("CCD finds renewed penetration while a rotating pair remains in contact",
          "[collision][ccd][angular][overlap]") {
    sqp::TestShape wall;
    wall.set_geometry(sqp::ConvexPoly({{-10, -0.5f}, {10, -0.5f}, {10, 0.5f}, {-10, 0.5f}}));
    const auto body = unit_box_shape();
    sqp::BodyState wall_start, wall_end, body_start, body_end;
    body_start.position = {0, sqp::real(1.5)};
    body_end = body_start;
    body_end.rotation = sqp::real(0.4);

    sqp::real fraction = 1;
    sqp::Vec2 normal;
    REQUIRE(sqp::sweep_shape_pair_overlap_approach_toi(wall, wall_start, wall_end, body, body_start,
                                                       body_end, fraction, normal));
    CHECK(fraction > 0);
    CHECK(fraction < 1);
    CHECK(finite_vec(normal));

    body_end.position.y = 3;
    CHECK_FALSE(sqp::sweep_shape_pair_overlap_approach_toi(wall, wall_start, wall_end, body,
                                                           body_start, body_end, fraction, normal));
    body_end = body_start;
    CHECK_FALSE(sqp::sweep_shape_pair_overlap_approach_toi(wall, wall_start, wall_end, body,
                                                           body_start, body_end, fraction, normal));
}

TEST_CASE("CCD overlap approach handles rotating offset circles",
          "[collision][ccd][angular][overlap]") {
    sqp::TestShape obstacle, moving;
    moving.set_geometry(sqp::Circle({sqp::real(0.5), 0}, 1));
    sqp::BodyState obstacle_start, obstacle_end, moving_start, moving_end;
    moving_start.position = {sqp::real(1.5), 0};
    moving_end = moving_start;
    moving_end.rotation = sqp::math::PI;

    for (bool polygon_obstacle : {false, true}) {
        if (polygon_obstacle)
            obstacle.set_geometry(sqp::ConvexPoly({{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}));
        else
            obstacle.set_geometry(sqp::Circle({}, 1));
        sqp::real fraction = 1;
        sqp::Vec2 normal;
        REQUIRE(sqp::sweep_shape_pair_overlap_approach_toi(obstacle, obstacle_start, obstacle_end,
                                                           moving, moving_start, moving_end,
                                                           fraction, normal));
        CHECK(fraction > 0);
        CHECK(fraction < 1);
        CHECK(finite_vec(normal));
        CHECK(normal.x < 0);
        REQUIRE(sqp::sweep_shape_pair_overlap_approach_toi(moving, moving_start, moving_end,
                                                           obstacle, obstacle_start, obstacle_end,
                                                           fraction, normal));
        CHECK(fraction > 0);
        CHECK(fraction < 1);
        CHECK(finite_vec(normal));
        CHECK(normal.x > 0);
    }
}

TEST_CASE("adaptive polygon sweep exposes a later approaching TOI",
          "[collision][ccd][shapecast][angular][regression]") {
    sqp::TestShape wall, body;
    wall.set_geometry(
        sqp::ConvexPoly({{-32.5f, -0.5f}, {32.5f, -0.5f}, {32.5f, 0.5f}, {-32.5f, 0.5f}}));
    body.set_geometry(sqp::ConvexPoly({{1.01400518f, 0.47668010f},
                                       {0.179995269f, 1.34693038f},
                                       {-0.834010005f, 0.870250046f},
                                       {-1.01400518f, -0.47668013f},
                                       {-0.179995105f, -1.34693050f},
                                       {0.834010184f, -0.870249808f}}));
    const sqp::Vec2 wall_position{130, 65};
    const sqp::Vec2 body_start{113.031105f, 63.0186195f};
    const sqp::Vec2 body_end{111.544907f, 64.8938293f};
    const sqp::real rotation_start = -12.7439938f;
    const sqp::real rotation_end = -13.4697332f;
    const auto first =
        sqp::convex_time_of_impact(wall, wall_position, 0, wall_position, 0, body, body_start,
                                   rotation_start, body_end, rotation_end, 0);
    REQUIRE(first.state == sqp::ConvexToiState::Hit);
    const auto first_position = body_start + (body_end - body_start) * first.fraction;
    const auto first_rotation = rotation_start + (rotation_end - rotation_start) * first.fraction;
    const auto kinematics = sqp::ccd_contact_kinematics(
        wall, wall_position, 0, {}, 0, body, first_position, first_rotation,
        {-380.46582f, 480.053925f}, -185.789368f, first.normal);
    CHECK(first.fraction == Catch::Approx(0.122297682f).margin(1e-5));
    CHECK(kinematics.normal_speed < 0);

    sqp::real toi = 1;
    sqp::Vec2 normal;
    REQUIRE(sqp::sweep_shape_pair_toi(wall, wall_position, 0, wall_position, 0, body, body_start,
                                      rotation_start, body_end, rotation_end, true, 0, toi,
                                      normal));
    const auto contact_position = body_start + (body_end - body_start) * toi;
    const auto contact_rotation = rotation_start + (rotation_end - rotation_start) * toi;
    const sqp::detail::ShapeQueryView wall_view = wall;
    const sqp::detail::ShapeQueryView body_view = body;
    const auto hit = sqp::sat_depth_query(wall_view.polygon_relative_to_body(), wall_position, 0,
                                          body_view.polygon_relative_to_body(), contact_position,
                                          contact_rotation);
    REQUIRE(hit.hit);
    CHECK(sqp::math::dot(normal, hit.normal) < sqp::real(-0.99));
}

TEST_CASE("rotating polygon finds wall reentry after a separating contact",
          "[collision][ccd][shapecast][reentry][regression]") {
    sqp::TestShape wall, body;
    wall.set_geometry(sqp::ConvexPoly({{-0.5f, -512}, {0.5f, -512}, {0.5f, 512}, {-0.5f, 512}}));
    body.set_geometry(sqp::ConvexPoly(
        {{2.47134066f, 2.35501337f}, {-4.23689795f, 0.276908308f}, {1.76555717f, -2.63192177f}}));
    const sqp::Vec2 wall_position{512.5f, 0};
    const sqp::Vec2 start{509.098f, -333.244f};
    const sqp::Vec2 end{511.358f, -328.538f};
    const sqp::real rotation_start = -29.0273f;
    const sqp::real rotation_end = -31.4326f;
    const sqp::detail::ShapeQueryView wall_view = wall;
    const sqp::detail::ShapeQueryView body_view = body;
    const auto wall_polygon = wall_view.polygon_relative_to_body();
    const auto body_polygon = body_view.polygon_relative_to_body();
    const auto overlap_at = [&](sqp::real fraction) {
        return sqp::sat_depth_query(wall_polygon, wall_position, 0, body_polygon,
                                    start + (end - start) * fraction,
                                    rotation_start + (rotation_end - rotation_start) * fraction)
            .hit;
    };
    REQUIRE(overlap_at(0));
    REQUIRE_FALSE(overlap_at(sqp::real(0.1)));
    REQUIRE(overlap_at(sqp::real(0.3)));
    sqp::real toi = 1;
    sqp::Vec2 normal;
    REQUIRE(sqp::sweep_shape_pair_reentry_toi(wall, wall_position, 0, wall_position, 0, body, start,
                                              rotation_start, end, rotation_end, true, 0, toi,
                                              normal));
    CHECK(toi > sqp::real(0.25));
    CHECK(toi < sqp::real(0.3));

    sqp::BodyState wall_start, wall_end, body_start, body_end;
    wall_start.position = wall_end.position = wall_position;
    body_start.position = start;
    body_start.rotation = rotation_start;
    body_end.position = end;
    body_end.rotation = rotation_end;
    REQUIRE(sqp::sweep_shape_pair_overlap_approach_toi(wall, wall_start, wall_end, body, body_start,
                                                       body_end, toi, normal));
    CHECK(toi > sqp::real(0.2));
    CHECK(toi < sqp::real(0.35));
}

TEST_CASE("convex CCD handles disabled settings invalid shapes and misses", "[collision][ccd]") {
    const sqp::TestShape a = unit_box_shape();
    sqp::TestShape b = unit_box_shape();
    sqp::Vec2 position_a0;
    sqp::real rotation_a0 = 0;
    sqp::Vec2 position_a1;
    sqp::real rotation_a1 = 0;
    position_a0 = {0, 0};
    position_a1 = {10, 0};
    sqp::Vec2 position_b0;
    sqp::real rotation_b0 = 0;
    sqp::Vec2 position_b1;
    sqp::real rotation_b1 = 0;
    position_b0 = position_b1 = {6, 4};
    bool enabled = true;
    sqp::real target_penetration_m = 1;
    target_penetration_m = 0;
    sqp::real toi = 1;
    sqp::Vec2 normal;

    CHECK_FALSE(sqp::sweep_shape_pair_toi(a, position_a0, rotation_a0, position_a1, rotation_a1, b,
                                          position_b0, rotation_b0, position_b1, rotation_b1,
                                          enabled, target_penetration_m, toi, normal));

    enabled = false;
    position_b0 = position_b1 = {6, 0};
    CHECK_FALSE(sqp::sweep_shape_pair_toi(a, position_a0, rotation_a0, position_a1, rotation_a1, b,
                                          position_b0, rotation_b0, position_b1, rotation_b1,
                                          enabled, target_penetration_m, toi, normal));

    enabled = true;
    b.set_geometry(sqp::ConvexPoly());
    CHECK_FALSE(sqp::sweep_shape_pair_toi(a, position_a0, rotation_a0, position_a1, rotation_a1, b,
                                          position_b0, rotation_b0, position_b1, rotation_b1,
                                          enabled, target_penetration_m, toi, normal));
}
