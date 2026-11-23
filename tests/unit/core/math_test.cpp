// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include <catch_amalgamated.hpp>

#include <cstdint>
#include <limits>
#include <type_traits>

#include "sqphysics2d/detail/math.h"

namespace {

template <class T> constexpr T tolerance() {
    return T(64) * std::numeric_limits<T>::epsilon();
}

template <class T> void check_near(sqp::Vec2T<T> actual, sqp::Vec2T<T> expected) {
    CHECK(actual.x == Catch::Approx(expected.x).margin(tolerance<T>()));
    CHECK(actual.y == Catch::Approx(expected.y).margin(tolerance<T>()));
}

} // namespace

TEMPLATE_TEST_CASE("Vec2: construction and arithmetic", "[math][vec2]", sqp::real32, sqp::real64) {
    using V = sqp::Vec2T<TestType>;

    constexpr V a{TestType(2), TestType(-3)};
    constexpr V b{TestType(-4), TestType(5)};

    STATIC_REQUIRE((V{} == V{0, 0}));
    STATIC_REQUIRE((V::zero() == V{0, 0}));
    STATIC_REQUIRE((V::one() == V{1, 1}));
    STATIC_REQUIRE((+a == a));
    STATIC_REQUIRE((-a == V{-2, 3}));
    STATIC_REQUIRE((a + b == V{-2, 2}));
    STATIC_REQUIRE((a - b == V{6, -8}));
    STATIC_REQUIRE((a * TestType(2) == V{4, -6}));
    STATIC_REQUIRE((TestType(2) * a == V{4, -6}));
    STATIC_REQUIRE((a / TestType(2) == V{1, TestType(-1.5)}));

    V value = a;
    V* const address = &value;
    CHECK(&(value += b) == address);
    CHECK((value == V{-2, 2}));
    CHECK(&(value -= b) == address);
    CHECK(value == a);
    CHECK(&(value *= TestType(2)) == address);
    CHECK((value == V{4, -6}));
    CHECK(&(value /= TestType(2)) == address);
    CHECK(value == a);
    CHECK_FALSE(value == b);
}

TEST_CASE("Vec2: explicit conversion between scalar types", "[math][vec2]") {
    constexpr sqp::Vec2f f{1.5f, -2.25f};
    constexpr sqp::Vec2d d{f};
    STATIC_REQUIRE((d == sqp::Vec2d{1.5, -2.25}));
    STATIC_REQUIRE((std::is_constructible_v<sqp::Vec2d, sqp::Vec2f>));
    STATIC_REQUIRE_FALSE((std::is_convertible_v<sqp::Vec2f, sqp::Vec2d>));
}

TEMPLATE_TEST_CASE("Vec2: length and vector products", "[math][vec2]", sqp::real32, sqp::real64) {
    using V = sqp::Vec2T<TestType>;
    const V v{TestType(3), TestType(4)};

    CHECK(V::zero().length_squared() == TestType(0));
    CHECK(V::zero().length() == TestType(0));
    CHECK(v.length_squared() == TestType(25));
    CHECK(v.length() == TestType(5));
    CHECK(sqp::math::dot(v, V{TestType(-2), TestType(5)}) == TestType(14));
    CHECK(sqp::math::dot(v, V{-4, 3}) == TestType(0));
    CHECK(sqp::math::cross(V{1, 0}, V{0, 1}) == TestType(1));
    CHECK(sqp::math::cross(V{0, 1}, V{1, 0}) == TestType(-1));
    CHECK(sqp::math::cross(v, v) == TestType(0));
    CHECK((sqp::math::cross(TestType(2), v) == V{-8, 6}));
    CHECK((sqp::math::cross(v, TestType(2)) == V{8, -6}));
    CHECK((sqp::math::perp_ccw(v) == V{-4, 3}));
}

TEMPLATE_TEST_CASE("Vec2: normalization and degenerate vectors", "[math][vec2]", sqp::real32,
                   sqp::real64) {
    using V = sqp::Vec2T<TestType>;
    const TestType eps = static_cast<TestType>(sqp::config::EPSILON);

    check_near(sqp::math::safe_normalized(V{3, 4}), V{TestType(0.6), TestType(0.8)});
    check_near(sqp::math::normalized_or_zero(V{-3, 4}), V{TestType(-0.6), TestType(0.8)});

    CHECK((sqp::math::safe_normalized(V{}) == V{1, 0}));
    CHECK((sqp::math::normalized_or_zero(V{}) == V{}));
    CHECK((sqp::math::safe_normalized(V{0, eps * TestType(0.5)}) == V{1, 0}));
    CHECK((sqp::math::normalized_or_zero(V{0, eps * TestType(0.5)}) == V{}));
    check_near(sqp::math::normalized_or_zero(V{0, eps * TestType(2)}), V{0, 1});
}

TEMPLATE_TEST_CASE("Vec2: rotation and edge normals", "[math][vec2]", sqp::real32, sqp::real64) {
    using V = sqp::Vec2T<TestType>;

    CHECK((sqp::math::rotate_vector(V{2, 3}, TestType(0), TestType(1)) == V{-3, 2}));
    CHECK((sqp::math::rotate_vector(V{2, 3}, TestType(1), TestType(0)) == V{2, 3}));
    check_near(sqp::math::rotate_vector(V{2, 0}, TestType(sqp::math::PI_REAL64 / 2)), V{0, 2});
    CHECK((sqp::math::outward_edge_normal_ccw(V{0, 0}, V{4, 0}) == V{0, -1}));
    CHECK((sqp::math::outward_edge_normal_ccw(V{0, 0}, V{0, 4}) == V{1, 0}));
    CHECK((sqp::math::outward_edge_normal_ccw(V{1, 1}, V{1, 1}) == V{1, 0}));
}

TEST_CASE("Aabb: empty, validity and inclusion", "[math][aabb]") {
    using namespace sqp;

    Aabb box = Aabb::empty();
    CHECK_FALSE(box.is_valid());
    box.include(Vec2{3, -2});
    CHECK(box.is_valid());
    CHECK((box == Aabb{Vec2{3, -2}, Vec2{3, -2}}));
    box.include(Vec2{-4, 5});
    CHECK((box == Aabb{Vec2{-4, -2}, Vec2{3, 5}}));
    box.include(Vec2{0, 0});
    CHECK((box == Aabb{Vec2{-4, -2}, Vec2{3, 5}}));

    CHECK((Aabb{Vec2{0, 0}, Vec2{0, 0}}.is_valid()));
    CHECK_FALSE((Aabb{Vec2{2, 0}, Vec2{1, 1}}.is_valid()));
    CHECK_FALSE((Aabb{Vec2{0, 2}, Vec2{1, 1}}.is_valid()));

    Aabb nonfinite{Vec2{0, 0}, Vec2{1, 1}};
    nonfinite.min.x = std::numeric_limits<real>::infinity();
    CHECK_FALSE(nonfinite.is_valid());
    nonfinite.min.x = std::numeric_limits<real>::quiet_NaN();
    CHECK_FALSE(nonfinite.is_valid());
    nonfinite.min.x = 0;
    nonfinite.max.y = -std::numeric_limits<real>::infinity();
    CHECK_FALSE(nonfinite.is_valid());
}

TEST_CASE("Aabb: geometry and expansion", "[math][aabb]") {
    using namespace sqp;
    Aabb box{Vec2{-2, 4}, Vec2{6, 10}};

    CHECK((box.center() == Vec2{2, 7}));
    CHECK((box.extent() == Vec2{4, 3}));
    box.expand(real(2));
    CHECK((box == Aabb{Vec2{-4, 2}, Vec2{8, 12}}));
    box.expand(real(0));
    CHECK((box == Aabb{Vec2{-4, 2}, Vec2{8, 12}}));
}

TEST_CASE("Aabb: inclusive overlap", "[math][aabb]") {
    using namespace sqp;
    const Aabb a{Vec2{0, 0}, Vec2{2, 2}};
    const Aabb intersect{Vec2{1, 1}, Vec2{3, 3}};
    const Aabb edge{Vec2{2, 0}, Vec2{4, 2}};
    const Aabb corner{Vec2{2, 2}, Vec2{4, 4}};
    const Aabb separate{Vec2{3, 0}, Vec2{4, 2}};
    const Aabb inside{Vec2{1, 1}, Vec2{1, 1}};

    CHECK(Aabb::overlap(a, a));
    CHECK(Aabb::overlap(a, intersect));
    CHECK(Aabb::overlap(intersect, a));
    CHECK(Aabb::overlap(a, edge));
    CHECK(Aabb::overlap(a, corner));
    CHECK(Aabb::overlap(a, inside));
    CHECK_FALSE(Aabb::overlap(a, separate));
    CHECK_FALSE(Aabb::overlap(separate, a));
}

TEST_CASE("Aabb: merge and swept bounds", "[math][aabb]") {
    using namespace sqp;
    const Aabb a{Vec2{1, 2}, Vec2{3, 4}};
    const Aabb b{Vec2{-2, 3}, Vec2{2, 8}};
    const Aabb merged{Vec2{-2, 2}, Vec2{3, 8}};

    CHECK((Aabb::merge(a, b) == merged));
    CHECK((Aabb::merge(b, a) == merged));
    CHECK((Aabb::merge(a, a) == a));
    CHECK((Aabb::merge(Aabb::empty(), a) == a));
    CHECK((Aabb::merge(a, Aabb::empty()) == a));

    CHECK((Aabb::swept(a, Vec2{0, 0}) == a));
    CHECK((Aabb::swept(a, Vec2{5, 6}) == Aabb{Vec2{1, 2}, Vec2{8, 10}}));
    CHECK((Aabb::swept(a, Vec2{-4, -6}) == Aabb{Vec2{-3, -4}, Vec2{3, 4}}));
    CHECK((Aabb::swept(a, Vec2{5, -6}) == Aabb{Vec2{1, -4}, Vec2{8, 4}}));
}

TEST_CASE("Math: angle conversions and wrapping", "[math][angle]") {
    using namespace sqp;

    CHECK(math::radian_to_degree(math::PI) == Catch::Approx(real(180)).margin(tolerance<real>()));
    CHECK(math::degree_to_radian(90) == Catch::Approx(math::PI / 2).margin(tolerance<real>()));
    CHECK(math::PI * math::INV_PI == Catch::Approx(real(1)).margin(tolerance<real>()));
    CHECK(math::DEG_PER_RAD * math::RAD_PER_DEG ==
          Catch::Approx(real(1)).margin(tolerance<real>()));
    CHECK(math::wrap_degree(real(180)) == real(-180));
    CHECK(math::wrap_degree(real(-180)) == real(-180));
    CHECK(math::wrap_degree(real(540)) == real(-180));
    CHECK(math::wrap_degree(real(-181)) == real(179));
    CHECK(math::wrap_radian(math::PI) == Catch::Approx(-math::PI).margin(tolerance<real>()));
    CHECK(math::wrap_radian(-math::PI) == Catch::Approx(-math::PI).margin(tolerance<real>()));
    CHECK(math::wrap(real(-1), real(0), real(5)) == real(4));
    CHECK(math::wrap(real(12), real(2), real(7)) == real(2));
    CHECK(math::wrap(real(9), real(4), real(4)) == real(4));
}

TEST_CASE("Math: scalar wrappers and clamp", "[math][scalar]") {
    using namespace sqp;

    CHECK(math::sqrt(real(9)) == real(3));
    CHECK(math::sin(real(0)) == real(0));
    CHECK(math::cos(real(0)) == real(1));
    CHECK(math::sin(math::PI / real(2)) == Catch::Approx(real(1)).margin(tolerance<real>()));
    CHECK(math::pow(real(2), 3) == real(8));
    CHECK(math::abs(real(-2.5)) == real(2.5));
    CHECK(math::floor(real(-1.2)) == real(-2));
    CHECK(math::ceil(real(-1.2)) == real(-1));
    CHECK(math::round(real(1.5)) == real(2));
    CHECK(math::round(real(-1.5)) == real(-2));
    CHECK(math::fmod(real(-5), real(2)) == real(-1));
    CHECK(math::isfinite(real(123)));
    CHECK_FALSE(math::isfinite(std::numeric_limits<real>::infinity()));
    CHECK_FALSE(math::isfinite(-std::numeric_limits<real>::infinity()));
    CHECK_FALSE(math::isfinite(std::numeric_limits<real>::quiet_NaN()));
    CHECK(math::clamp01(real(-0.5)) == real(0));
    CHECK(math::clamp01(real(0)) == real(0));
    CHECK(math::clamp01(real(0.25)) == real(0.25));
    CHECK(math::clamp01(real(1)) == real(1));
    CHECK(math::clamp01(real(1.5)) == real(1));
}

TEST_CASE("Math: time-unit conversions", "[math][time]") {
    using namespace sqp;
    constexpr real units = real(config::TU_PER_SECOND);

    STATIC_REQUIRE(math::tu_to_seconds<real>(config::TU_PER_SECOND) == real(1));
    STATIC_REQUIRE(math::seconds_to_tu(real(-1)) == tick_t(0));
    STATIC_REQUIRE(math::seconds_to_tu(real(0)) == tick_t(0));
    STATIC_REQUIRE(math::seconds_to_tu(real(1)) == tick_t(config::TU_PER_SECOND));
    CHECK(math::seconds_to_tu(real(2.25) / units) == tick_t(2));
    CHECK(math::seconds_to_tu(real(2.75) / units) == tick_t(3));
    CHECK(math::tu_to_seconds<real>(tick_t(0)) == real(0));
    CHECK(math::tu_to_seconds<real64>(config::TU_PER_SECOND) == 1.0);
}

TEST_CASE("Math: bit flags", "[math][flags]") {
    using namespace sqp;
    enum class Flag : std::uint32_t { A = 1, B = 2, C = 4 };

    std::uint32_t flags = 0;
    CHECK_FALSE(math::has_any_flags(flags, Flag::A));
    CHECK_FALSE(math::has_all_flags(flags, Flag::A));
    CHECK_FALSE(math::has_any_flags(flags, 0u));
    CHECK(math::has_all_flags(flags, 0u));

    math::enable_flags(flags, Flag::A);
    math::enable_flags(flags, Flag::C);
    CHECK(flags == 5u);
    CHECK(math::has_any_flags(flags, 3u));
    CHECK(math::has_all_flags(flags, Flag::A));
    CHECK(math::has_all_flags(flags, 5u));
    CHECK_FALSE(math::has_all_flags(flags, 3u));

    math::disable_flags(flags, Flag::A);
    CHECK(flags == 4u);
    math::set_flags(flags, Flag::B, true);
    CHECK(flags == 6u);
    math::set_flags(flags, Flag::B, false);
    CHECK(flags == 4u);
    math::set_flags(flags, Flag::C, true);
    CHECK(flags == 4u);
}

TEST_CASE("Math: floor_div with negative coordinates", "[math][grid]") {
    using namespace sqp;

    CHECK(math::floor_div(real(0), real(4)) == 0);
    CHECK(math::floor_div(real(3.9), real(4)) == 0);
    CHECK(math::floor_div(real(4), real(4)) == 1);
    CHECK(math::floor_div(real(-0.1), real(4)) == -1);
    CHECK(math::floor_div(real(-4), real(4)) == -1);
    CHECK(math::floor_div(real(-4.1), real(4)) == -2);
    CHECK(math::floor_div(real(9), real(4)) == 2);
}
