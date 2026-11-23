#include "ccd/convex_toi.h"
#include "ccd/shape_ccd.h"
#include "sqphysics2d/detail/data_storage.h"

#include <catch_amalgamated.hpp>

namespace {

sqp::PolygonShapeDef box(sqp::Vec2 half, sqp::Vec2 offset = {}) {
    sqp::PolygonShapeDef def;
    def.vertices = {{-half.x, -half.y}, {half.x, -half.y}, {half.x, half.y}, {-half.x, half.y}};
    def.local_offset = offset;
    return def;
}

} // namespace

TEST_CASE("shape pool queries circle and polygon geometry", "[shape][query]") {
    sqp::ShapePool pool;
    const auto circle = pool.add_circle({2});
    const auto polygon = pool.add_polygon(box({3, 4}));
    REQUIRE(pool.is_valid(circle));
    REQUIRE(pool.is_valid(polygon));
    CHECK(pool.sweep_radius(circle) == Catch::Approx(2));
    CHECK(pool.sweep_radius(polygon) == Catch::Approx(5));

    sqp::BodyState state;
    state.position = {3, 4};
    const auto bounds = pool.world_aabb(polygon, state);
    CHECK(bounds.min.x == Catch::Approx(0));
    CHECK(bounds.min.y == Catch::Approx(0));
    CHECK(bounds.max.x == Catch::Approx(6));
    CHECK(bounds.max.y == Catch::Approx(8));
}

TEST_CASE("shape pool swept bounds contain rotated geometry", "[shape][query][aabb]") {
    sqp::ShapePool pool;
    const auto shape = pool.add_polygon(box({1, 2}));
    sqp::BodyState start;
    start.position = {3, 4};
    auto end = start;
    end.rotation = sqp::math::PI / 2;
    const auto swept = pool.swept_world_aabb(shape, start, end);
    CHECK(swept.min.x == Catch::Approx(3 - sqp::math::sqrt(sqp::real(5))).margin(1e-5));
    CHECK(swept.max.y == Catch::Approx(4 + sqp::math::sqrt(sqp::real(5))).margin(1e-5));
    for (int i = 0; i <= 16; ++i) {
        auto sample = start;
        sample.rotation = end.rotation * sqp::real(i) / 16;
        const auto bounds = pool.world_aabb(shape, sample);
        CHECK(bounds.min.x >= swept.min.x - sqp::config::EPSILON);
        CHECK(bounds.min.y >= swept.min.y - sqp::config::EPSILON);
        CHECK(bounds.max.x <= swept.max.x + sqp::config::EPSILON);
        CHECK(bounds.max.y <= swept.max.y + sqp::config::EPSILON);
    }
}

TEST_CASE("shape pool offsets affect bounds and ray casts", "[shape][query][raycast]") {
    sqp::ShapePool pool;
    const auto circle = pool.add_circle({1, {3, 0}});
    const auto polygon = pool.add_polygon(box({1, 1}, {3, 0}));
    sqp::BodyState state;
    state.position = {10, 20};
    state.rotation = sqp::math::PI / 2;
    for (auto shape : {circle, polygon}) {
        const auto bounds = pool.world_aabb(shape, state);
        CHECK(bounds.min.x == Catch::Approx(9).margin(1e-5));
        CHECK(bounds.max.x == Catch::Approx(11).margin(1e-5));
        CHECK(bounds.min.y == Catch::Approx(22).margin(1e-5));
        CHECK(bounds.max.y == Catch::Approx(24).margin(1e-5));
        sqp::real fraction = 1;
        sqp::Vec2 normal;
        REQUIRE(pool.ray_cast(shape, state, {6, 23}, {8, 0}, fraction, normal));
        CHECK(fraction == Catch::Approx(sqp::real(0.375)).margin(1e-5));
        CHECK(normal.x == Catch::Approx(-1).margin(1e-5));
    }
}

TEST_CASE("shape pool queries use updated geometry and reject stale handles", "[shape][query]") {
    sqp::ShapePool pool;
    const auto old = pool.add_circle({1, {2, 0}});
    const auto polygon = pool.add_polygon(box({1, 1}));
    sqp::BodyState state;
    CHECK(pool.sweep_radius(old) == Catch::Approx(3));
    REQUIRE(pool.update_circle(old, {2, {4, 0}}));
    CHECK(pool.sweep_radius(old) == Catch::Approx(6));
    CHECK(pool.world_aabb(old, state).min.x == Catch::Approx(2));
    REQUIRE(pool.update_polygon(polygon, box({2, 1}, {5, 0})));
    CHECK(pool.world_aabb(polygon, state).min.x == Catch::Approx(3));
    REQUIRE(pool.remove(old));
    const auto replacement = pool.add_circle({1});
    REQUIRE(pool.is_valid(replacement));
    CHECK_FALSE(pool.is_valid(old));
    CHECK_FALSE(pool.world_aabb(old, state).is_valid());
    CHECK_FALSE(pool.swept_world_aabb(old, state, state).is_valid());
    CHECK(pool.sweep_radius(old) == 0);
    sqp::real fraction = 1;
    sqp::Vec2 normal;
    CHECK_FALSE(pool.ray_cast(old, state, {}, {10, 0}, fraction, normal));
}

TEST_CASE("CCD pair queries read handles from one shape pool", "[shape][query][ccd]") {
    sqp::ShapePool pool;
    const auto a = pool.add_circle({1});
    const auto b = pool.add_circle({1});
    const auto distance = sqp::convex_signed_distance(pool, a, {0, 0}, 0, b, {4, 0}, 0);
    REQUIRE(distance.valid);
    CHECK(distance.separation == Catch::Approx(2));
    sqp::real toi = 1;
    sqp::Vec2 normal;
    REQUIRE(sqp::sweep_shape_pair_toi(pool, a, {0, 0}, 0, {4, 0}, 0, b, {4, 0}, 0, {4, 0}, 0, true,
                                      0, toi, normal));
    CHECK(toi == Catch::Approx(sqp::real(0.5)).margin(1e-5));
    REQUIRE(pool.remove(b));
    CHECK_FALSE(sqp::sweep_shape_pair_toi(pool, a, {0, 0}, 0, {4, 0}, 0, b, {4, 0}, 0, {4, 0}, 0,
                                          true, 0, toi, normal));
}
