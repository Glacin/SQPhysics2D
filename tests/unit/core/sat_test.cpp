#include "narrow_phase/sat.h"
#include <catch_amalgamated.hpp>
#include <random>

namespace {

sqp::ConvexPoly box(sqp::real half_width = 1, sqp::real half_height = 1) {
    return sqp::ConvexPoly({{-half_width, -half_height},
                            {half_width, -half_height},
                            {half_width, half_height},
                            {-half_width, half_height}});
}

void check_sat_queries_agree(const sqp::ConvexPoly& a, sqp::Vec2 position_a, sqp::real rotation_a,
                             const sqp::ConvexPoly& b, sqp::Vec2 position_b, sqp::real rotation_b,
                             bool expected_hit) {
    const auto world_a = sqp::transform_poly_to_world(a, position_a, rotation_a);
    const auto world_b = sqp::transform_poly_to_world(b, position_b, rotation_b);
    const auto full =
        sqp::sat_intersect(a, position_a, rotation_a, b, position_b, rotation_b, world_a, world_b);
    const auto depth = sqp::sat_depth_query(a, position_a, rotation_a, b, position_b, rotation_b);

    REQUIRE(full.hit == expected_hit);
    CHECK(depth.hit == full.hit);
    CHECK(depth.ref_edge_index == -1);
    if (!full.hit)
        return;

    CHECK(depth.depth == Catch::Approx(full.depth).margin(1e-5));
    CHECK(depth.normal.x == Catch::Approx(full.normal.x).margin(1e-5));
    CHECK(depth.normal.y == Catch::Approx(full.normal.y).margin(1e-5));
    CHECK(depth.ref_is_A == full.ref_is_A);
    CHECK(full.ref_edge_index >= 0);
    CHECK(full.ref_edge_index <
          static_cast<sqp::i32>(full.ref_is_A ? world_a.size() : world_b.size()));
}

} // namespace

TEST_CASE("SAT reports overlap depth and A to B normal", "[narrow_phase][sat]") {
    const sqp::SatHit hit = sqp::sat_depth_query(box(), {0, 0}, 0, box(), {sqp::real(1.5), 0}, 0);

    REQUIRE(hit.hit);
    CHECK(hit.depth == Catch::Approx(0.5));
    CHECK(hit.normal.x == Catch::Approx(1));
    CHECK(hit.normal.y == Catch::Approx(0));
}

TEST_CASE("SAT rejects separated and invalid polygons", "[narrow_phase][sat]") {
    CHECK_FALSE(sqp::sat_depth_query(box(), {0, 0}, 0, box(), {3, 0}, 0).hit);
    CHECK_FALSE(sqp::sat_depth_query(sqp::ConvexPoly(), {0, 0}, 0, box(), {0, 0}, 0).hit);
}

TEST_CASE("SAT treats touching edges as zero-depth contact", "[narrow_phase][sat]") {
    const sqp::SatHit hit = sqp::sat_depth_query(box(), {0, 0}, 0, box(), {2, 0}, 0);

    REQUIRE(hit.hit);
    CHECK(hit.depth == Catch::Approx(0).margin(1e-6));
    CHECK(hit.normal.x == Catch::Approx(1));
}

TEST_CASE("SAT containment depth selects the nearest wall face", "[narrow_phase][sat]") {
    const auto wall = box(129, sqp::real(0.5));
    const auto body = box(sqp::real(0.15), sqp::real(0.153));
    const sqp::SatHit hit = sqp::sat_depth_query(wall, {0, sqp::real(-128.5)}, 0, body,
                                                 {sqp::real(127.849), sqp::real(-127.848)}, 0);

    REQUIRE(hit.hit);
    CHECK(hit.depth == Catch::Approx(0.001).margin(2e-4));
    CHECK(hit.normal.x == Catch::Approx(0).margin(1e-6));
    CHECK(hit.normal.y == Catch::Approx(1));
}

TEST_CASE("SAT handles rotated boxes with stable normal direction", "[narrow_phase][sat]") {
    const sqp::Vec2 offset{sqp::real(1.2), sqp::real(0.15)};
    const sqp::SatHit hit = sqp::sat_depth_query(box(), {0, 0}, sqp::math::degree_to_radian(15),
                                                 box(), offset, sqp::math::degree_to_radian(-20));

    REQUIRE(hit.hit);
    CHECK(hit.depth > sqp::real(0));
    CHECK(sqp::math::dot(hit.normal, offset) > sqp::real(0));
}

TEST_CASE("depth-only SAT agrees with the reference SAT on rotated boxes",
          "[narrow_phase][sat][ccd]") {
    std::mt19937 rng(0x51a7u);
    std::uniform_real_distribution<float> coordinate(-3.0f, 3.0f);
    std::uniform_real_distribution<float> angle(-3.0f, 3.0f);
    const auto a = box(sqp::real(1.25), sqp::real(0.4));
    const auto b = box(sqp::real(0.65), sqp::real(1.1));
    for (int sample = 0; sample < 1000; ++sample) {
        const sqp::Vec2 offset{coordinate(rng), coordinate(rng)};
        const sqp::real rotation_a = angle(rng), rotation_b = angle(rng);
        const auto world_a = sqp::transform_poly_to_world(a, {}, rotation_a);
        const auto world_b = sqp::transform_poly_to_world(b, offset, rotation_b);
        const auto reference =
            sqp::sat_intersect(a, {}, rotation_a, b, offset, rotation_b, world_a, world_b);
        const auto depth = sqp::sat_depth_query(a, {}, rotation_a, b, offset, rotation_b);
        REQUIRE(depth.hit == reference.hit);
        if (reference.hit) {
            CHECK(depth.depth == Catch::Approx(reference.depth).margin(1e-5));
            CHECK(depth.normal.x == Catch::Approx(reference.normal.x).margin(1e-5));
            CHECK(depth.normal.y == Catch::Approx(reference.normal.y).margin(1e-5));
        }
    }
}

TEST_CASE("SAT query variants agree for touching, separation, containment and offset polygons",
          "[narrow_phase][sat]") {
    const auto unit = box();
    check_sat_queries_agree(unit, {}, 0, unit, {2, 0}, 0, true);
    check_sat_queries_agree(unit, {}, 0, unit, {3, 0}, 0, false);
    check_sat_queries_agree(box(4, 2), {}, 0, box(sqp::real(0.5)), {sqp::real(3.75), 0}, 0, true);
    check_sat_queries_agree(unit, {}, sqp::math::degree_to_radian(15), unit,
                            {sqp::real(1.2), sqp::real(0.15)}, sqp::math::degree_to_radian(-20),
                            true);
    const sqp::ConvexPoly offset({{2, 0}, {4, 0}, {3, 2}});
    check_sat_queries_agree(offset, {}, sqp::real(0.35), offset, {sqp::real(0.4), sqp::real(0.2)},
                            sqp::real(-0.1), true);
}

TEST_CASE("polygon transform applies body position and centroid offset", "[narrow_phase][sat]") {
    const sqp::ConvexPoly poly({{2, 0}, {4, 0}, {3, 2}});
    const sqp::PolygonVertices world = sqp::transform_poly_to_world(poly, {10, 0}, 0);
    const sqp::PolygonVertices rotated =
        sqp::transform_poly_to_world(poly, {10, 0}, sqp::math::PI / 2);

    REQUIRE(world.size() == 3);
    REQUIRE(rotated.size() == 3);
    CHECK(world[0].x == Catch::Approx(12));
    CHECK(rotated[0].x == Catch::Approx(10).margin(1e-5));
    CHECK(rotated[0].y == Catch::Approx(2).margin(1e-5));
}

TEST_CASE("polygon validation uses fixed-capacity vertices", "[narrow_phase][sat][polygon]") {
    const sqp::ConvexPoly cleaned({{0, 0}, {1, 0}, {2, 0}, {2, 2}, {0, 2}, {0, 0}});
    REQUIRE(cleaned.local.size() >= 3);
    CHECK(cleaned.vertex_count() == 4);

    sqp::PolygonVertices oversized;
    for (std::size_t i = 0; i < sqp::config::MAX_POLYGON_VERTICES + 1; ++i) {
        oversized.push_back({static_cast<sqp::real>(i), 0});
    }
    CHECK(oversized.overflowed());
    CHECK(sqp::ConvexPoly(oversized).local.size() < 3);
}
