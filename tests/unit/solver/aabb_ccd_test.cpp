#include "catch_amalgamated.hpp"
#include "ccd/shape_ccd.h"

TEST_CASE("circle versus AABB CCD uses relative translation", "[collision][ccd]") {
    bool enabled = true;
    sqp::real target_penetration_m = 0;
    sqp::real toi = 1;
    sqp::Vec2 normal;
    REQUIRE(sqp::sweep_circle_vs_aabb_toi({0, 0}, {10, 0}, 1, {{5, -1}, {7, 1}}, {{7, -1}, {9, 1}},
                                          enabled, target_penetration_m, toi, normal));
    CHECK(toi == Catch::Approx(0.5f));
    CHECK(normal.x == Catch::Approx(-1));
    CHECK(normal.y == Catch::Approx(0));
}

TEST_CASE("circle versus AABB CCD handles initial overlap and tangent miss", "[collision][ccd]") {
    bool enabled = true;
    sqp::real target_penetration_m = 0;
    sqp::real toi = 1;
    sqp::Vec2 normal;
    REQUIRE(sqp::sweep_circle_vs_aabb_toi({0, 1.5f}, {1, 1.5f}, 0.5f, {{-1, -1}, {1, 1}},
                                          {{-1, -1}, {1, 1}}, enabled, target_penetration_m, toi,
                                          normal));
    CHECK(toi == Catch::Approx(0));

    CHECK_FALSE(sqp::sweep_circle_vs_aabb_toi({-3, 1.5001f}, {3, 1.5001f}, 0.5f, {{-1, -1}, {1, 1}},
                                              {{-1, -1}, {1, 1}}, enabled, target_penetration_m,
                                              toi, normal));
}

TEST_CASE("circle versus AABB CCD rejects disabled invalid bounds and zero-radius cases",
          "[collision][ccd]") {
    bool enabled = true;
    sqp::real target_penetration_m = 0;
    sqp::real toi = 1;
    sqp::Vec2 normal;

    enabled = false;
    CHECK_FALSE(sqp::sweep_circle_vs_aabb_toi({0, 0}, {10, 0}, 1, {{5, -1}, {7, 1}},
                                              {{5, -1}, {7, 1}}, enabled, target_penetration_m, toi,
                                              normal));

    enabled = true;
    CHECK_FALSE(sqp::sweep_circle_vs_aabb_toi({0, 0}, {10, 0}, 0, {{5, -1}, {7, 1}},
                                              {{5, -1}, {7, 1}}, enabled, target_penetration_m, toi,
                                              normal));
    CHECK_FALSE(sqp::sweep_circle_vs_aabb_toi({0, 0}, {10, 0}, 1, {{7, 1}, {5, -1}},
                                              {{7, 1}, {5, -1}}, enabled, target_penetration_m, toi,
                                              normal));
}

TEST_CASE("circle versus AABB CCD tolerates float extent rounding during translation",
          "[collision][ccd]") {
    sqp::real toi = 1;
    sqp::Vec2 normal;
    const sqp::Aabb box0{{126.377922058f, 64.045356750f},
                         {126.847221375f, 64.902763367f}};
    const sqp::Aabb box1{{130.662094116f, 63.094112396f},
                         {131.131393433f, 63.951511383f}};
    REQUIRE(sqp::sweep_circle_vs_aabb_toi({128.5f, 64}, {128.5f, 64}, 0.2f,
                                          box0, box1, true, 0, toi, normal));
    CHECK(toi > 0);
    CHECK(toi < 1);
}

TEST_CASE("circle versus AABB CCD reaches rounded corner at the correct time", "[collision][ccd]") {
    sqp::real toi;
    sqp::Vec2 normal;
    const sqp::Aabb box{{-1, -1}, {1, 1}};
    REQUIRE(sqp::sweep_circle_vs_aabb_toi({-3, sqp::real(1.4)}, {0, sqp::real(1.4)}, sqp::real(0.5),
                                          box, box, true, 0, toi, normal));
    CHECK(toi == Catch::Approx(sqp::real(1.7) / sqp::real(3)).margin(1e-5));
    CHECK(normal.x == Catch::Approx(sqp::real(-0.6)).margin(1e-4));
    CHECK(normal.y == Catch::Approx(sqp::real(0.8)).margin(1e-4));
}

TEST_CASE("AABB versus AABB CCD reports face TOI", "[collision][ccd]") {
    bool enabled = true;
    sqp::real target_penetration_m = 0;
    sqp::real toi = 1;
    sqp::Vec2 normal;
    REQUIRE(sqp::sweep_aabb_vs_aabb_toi({{-1, -1}, {1, 1}}, {{9, -1}, {11, 1}}, {{5, -2}, {7, 2}},
                                        {{5, -2}, {7, 2}}, enabled, target_penetration_m, toi,
                                        normal));
    CHECK(toi == Catch::Approx(0.4f));
    CHECK(normal.x == Catch::Approx(-1));
    CHECK(normal.y == Catch::Approx(0));
}

TEST_CASE("AABB versus AABB CCD handles relative motion and misses", "[collision][ccd]") {
    bool enabled = true;
    sqp::real target_penetration_m = 0;
    sqp::real toi = 1;
    sqp::Vec2 normal;
    REQUIRE(sqp::sweep_aabb_vs_aabb_toi({{-1, -1}, {1, 1}}, {{5, -1}, {7, 1}}, {{7, -1}, {9, 1}},
                                        {{5, -1}, {7, 1}}, enabled, target_penetration_m, toi,
                                        normal));
    CHECK(toi == Catch::Approx(0.75f));
    CHECK(normal.x == Catch::Approx(-1));

    CHECK_FALSE(sqp::sweep_aabb_vs_aabb_toi({{-1, -1}, {1, 1}}, {{5, -1}, {7, 1}}, {{7, 3}, {9, 5}},
                                            {{5, 3}, {7, 5}}, enabled, target_penetration_m, toi,
                                            normal));
}

TEST_CASE("AABB versus AABB CCD respects target penetration and invalid bounds",
          "[collision][ccd]") {
    bool enabled = true;
    sqp::real target_penetration_m = 0;
    target_penetration_m = sqp::real(0.5);
    sqp::real toi = 1;
    sqp::Vec2 normal;
    REQUIRE(sqp::sweep_aabb_vs_aabb_toi({{-1, -1}, {1, 1}}, {{9, -1}, {11, 1}}, {{5, -2}, {7, 2}},
                                        {{5, -2}, {7, 2}}, enabled, target_penetration_m, toi,
                                        normal));
    CHECK(toi == Catch::Approx(0.45f));
    CHECK(normal.x == Catch::Approx(-1));

    CHECK_FALSE(sqp::sweep_aabb_vs_aabb_toi({{1, 1}, {-1, -1}}, {{9, -1}, {11, 1}},
                                            {{5, -2}, {7, 2}}, {{5, -2}, {7, 2}}, enabled,
                                            target_penetration_m, toi, normal));
}

TEST_CASE("AABB versus AABB CCD tolerates float extent rounding during translation",
          "[collision][ccd]") {
    sqp::real toi = 1;
    sqp::Vec2 normal;
    const sqp::Aabb box0{{126.377922058f, 64.045356750f},
                         {126.847221375f, 64.902763367f}};
    const sqp::Aabb box1{{130.662094116f, 63.094112396f},
                         {131.131393433f, 63.951511383f}};
    const sqp::Aabb wall{{128, -128}, {129, 128}};
    REQUIRE(sqp::sweep_aabb_vs_aabb_toi(box0, box1, wall, wall, true, 0, toi, normal));
    CHECK(toi == Catch::Approx(0.2690785f).margin(1e-4f));
    CHECK(normal.x == Catch::Approx(-1));
    CHECK(normal.y == Catch::Approx(0));
}
