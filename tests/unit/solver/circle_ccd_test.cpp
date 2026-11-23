#include "catch_amalgamated.hpp"
#include "ccd/shape_ccd.h"

namespace {
sqp::ConvexPoly unit_box() {
    return sqp::ConvexPoly({{-1, -1}, {1, -1}, {1, 1}, {-1, 1}});
}
} // namespace

TEST_CASE("circle CCD reports contact", "[collision][ccd]") {
    bool enabled = true;
    sqp::real target_penetration_m = 1;
    target_penetration_m = 0;
    sqp::real toi = 1;
    sqp::Vec2 normal;
    REQUIRE(sqp::sweep_circle_vs_circle_toi({0, 0}, {10, 0}, 1, {6, 0}, {6, 0}, 1, enabled,
                                            target_penetration_m, toi, normal));
    CHECK(toi == Catch::Approx(0.4f));
    CHECK(normal.x == Catch::Approx(-1));
}

TEST_CASE("circle CCD handles relative motion initial overlap and misses", "[collision][ccd]") {
    bool enabled = true;
    sqp::real target_penetration_m = 1;
    target_penetration_m = 0;
    sqp::real toi = 1;
    sqp::Vec2 normal;

    REQUIRE(sqp::sweep_circle_vs_circle_toi({0, 0}, {6, 0}, 1, {8, 0}, {4, 0}, 1, enabled,
                                            target_penetration_m, toi, normal));
    CHECK(toi == Catch::Approx(sqp::real(0.6)).margin(1e-5));
    CHECK(normal.x == Catch::Approx(-1));

    REQUIRE(sqp::sweep_circle_vs_circle_toi({0, 0}, {1, 0}, 2, {1, 0}, {2, 0}, 2, enabled,
                                            target_penetration_m, toi, normal));
    CHECK(toi == Catch::Approx(0));

    CHECK_FALSE(sqp::sweep_circle_vs_circle_toi({0, 0}, {0, 0}, 1, {5, 0}, {5, 0}, 1, enabled,
                                                target_penetration_m, toi, normal));
    CHECK_FALSE(sqp::sweep_circle_vs_circle_toi({0, 2.0001f}, {10, 2.0001f}, 1, {5, 0}, {5, 0}, 1,
                                                enabled, target_penetration_m, toi, normal));
}

TEST_CASE("circle CCD respects disabled settings target penetration and zero motion",
          "[collision][ccd]") {
    bool enabled = true;
    sqp::real target_penetration_m = 1;
    target_penetration_m = 0;
    sqp::real toi = 1;
    sqp::Vec2 normal;

    enabled = false;
    CHECK_FALSE(sqp::sweep_circle_vs_circle_toi({0, 0}, {10, 0}, 1, {6, 0}, {6, 0}, 1, enabled,
                                                target_penetration_m, toi, normal));

    enabled = true;
    CHECK_FALSE(sqp::sweep_circle_vs_circle_toi({0, 0}, {0, 0}, 1, {6, 0}, {6, 0}, 1, enabled,
                                                target_penetration_m, toi, normal));
    CHECK_FALSE(sqp::sweep_circle_vs_circle_toi({0, 0}, {10, 0}, 0, {6, 0}, {6, 0}, 0, enabled,
                                                target_penetration_m, toi, normal));

    target_penetration_m = sqp::real(0.5);
    REQUIRE(sqp::sweep_circle_vs_circle_toi({0, 0}, {10, 0}, 1, {6, 0}, {6, 0}, 1, enabled,
                                            target_penetration_m, toi, normal));
    CHECK(toi == Catch::Approx(0.45f));
    CHECK(normal.x == Catch::Approx(-1));
}

TEST_CASE("circle CCD hits polygon face", "[collision][ccd]") {
    const sqp::ConvexPoly box({{-5, -5}, {5, -5}, {5, 5}, {-5, 5}});
    bool enabled = true;
    sqp::real target_penetration_m = 1;
    target_penetration_m = 0.05f;
    sqp::real toi = 1;
    sqp::Vec2 normal;
    sqp::Vec2 point;
    REQUIRE(sqp::sweep_circle_vs_poly_toi({19.1f, 2.5f}, {3.475f, 4.0625f}, 0.25f, box, {}, 0, 0,
                                          {}, enabled, target_penetration_m, toi, normal, point));
    CHECK(toi == Catch::Approx(0.8896f).margin(1e-5));
    CHECK(normal.x == Catch::Approx(1));
    CHECK(point.x == Catch::Approx(5));
}

TEST_CASE("circle CCD polygon face hits are angle stable") {
    bool enabled = true;
    sqp::real target_penetration_m = 1;
    target_penetration_m = 0;
    bool all_face_sweeps_match = true;
    sqp::real first_failed_y = 0;
    for (const sqp::real y : {-0.9f, -0.45f, 0.0f, 0.45f, 0.9f}) {
        sqp::real toi = 1;
        sqp::Vec2 normal;
        sqp::Vec2 point;
        const bool hit =
            sqp::sweep_circle_vs_poly_toi({-10, y}, {10, y}, 0.5f, unit_box(), {}, 0, 0, {},
                                          enabled, target_penetration_m, toi, normal, point);
        const bool matches = hit && toi == Catch::Approx(0.425f).margin(1e-5) &&
                             normal.x == Catch::Approx(-1) && normal.y == Catch::Approx(0);
        if (!matches && all_face_sweeps_match)
            first_failed_y = y;
        all_face_sweeps_match = matches && all_face_sweeps_match;
    }
    CAPTURE(first_failed_y);
    CHECK(all_face_sweeps_match);
}

TEST_CASE("circle CCD polygon corner hits are angle stable") {
    bool enabled = true;
    sqp::real target_penetration_m = 1;
    target_penetration_m = 0;
    const sqp::Vec2 corner(1, 1);
    bool all_corner_sweeps_match = true;
    sqp::real first_failed_degrees = 0;
    for (const sqp::real degrees : {5.0f, 20.0f, 35.0f, 45.0f, 55.0f, 70.0f, 85.0f}) {
        const sqp::real angle = degrees * sqp::math::PI / 180;
        const sqp::Vec2 direction(sqp::math::cos(angle), sqp::math::sin(angle));
        sqp::real toi = 1;
        sqp::Vec2 normal;
        sqp::Vec2 point;
        const bool hit = sqp::sweep_circle_vs_poly_toi(
            corner + direction * 20, corner - direction * 5, 0.5f, unit_box(), {}, 0, 0, {},
            enabled, target_penetration_m, toi, normal, point);
        const bool matches = hit && toi == Catch::Approx(0.78f).margin(1e-5) &&
                             normal.x == Catch::Approx(direction.x).margin(1e-5) &&
                             normal.y == Catch::Approx(direction.y).margin(1e-5) &&
                             point.x == Catch::Approx(1).margin(1e-5) &&
                             point.y == Catch::Approx(1).margin(1e-5);
        if (!matches && all_corner_sweeps_match)
            first_failed_degrees = degrees;
        all_corner_sweeps_match = matches && all_corner_sweeps_match;
    }
    CAPTURE(first_failed_degrees);
    CHECK(all_corner_sweeps_match);
}

TEST_CASE("circle CCD polygon near tangent misses") {
    bool enabled = true;
    sqp::real target_penetration_m = 1;
    target_penetration_m = 0;
    sqp::real toi = 1;
    sqp::Vec2 normal;
    sqp::Vec2 point;
    CHECK_FALSE(sqp::sweep_circle_vs_poly_toi({-10, 1.5001f}, {10, 1.5001f}, 0.5f, unit_box(), {},
                                              0, 0, {}, enabled, target_penetration_m, toi, normal,
                                              point));
    CHECK_FALSE(sqp::sweep_circle_vs_poly_toi({1.49f, 1.49f}, {3, 3}, 0.5f, unit_box(), {}, 0, 0,
                                              {}, enabled, target_penetration_m, toi, normal,
                                              point));
}

TEST_CASE("circle polygon CCD handles invalid input disabled settings and initial overlap normal",
          "[collision][ccd]") {
    bool enabled = true;
    sqp::real target_penetration_m = 1;
    target_penetration_m = 0;
    sqp::real toi = 1;
    sqp::Vec2 normal;
    sqp::Vec2 point;

    enabled = false;
    CHECK_FALSE(sqp::sweep_circle_vs_poly_toi({-10, 0}, {10, 0}, 0.5f, unit_box(), {}, 0, 0, {},
                                              enabled, target_penetration_m, toi, normal, point));

    enabled = true;
    CHECK_FALSE(sqp::sweep_circle_vs_poly_toi({-10, 0}, {10, 0}, 0, unit_box(), {}, 0, 0, {},
                                              enabled, target_penetration_m, toi, normal, point));
    CHECK_FALSE(sqp::sweep_circle_vs_poly_toi({-10, 0}, {10, 0}, 0.5f, sqp::ConvexPoly(), {}, 0, 0,
                                              {}, enabled, target_penetration_m, toi, normal,
                                              point));

    REQUIRE(sqp::sweep_circle_vs_poly_toi({0, 0}, {0, 0}, 0.5f, unit_box(), {}, 0, 0, {}, enabled,
                                          target_penetration_m, toi, normal, point));
    CHECK(toi == Catch::Approx(0));
    CHECK(sqp::math::isfinite(normal.x));
    CHECK(sqp::math::isfinite(normal.y));
}
