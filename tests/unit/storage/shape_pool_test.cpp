#include <array>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

#include "sqphysics2d/detail/data_storage.h"
#include <catch_amalgamated.hpp>

namespace {

sqp::PolygonShapeDef square(sqp::real half_extent) {
    return {{{-half_extent, -half_extent},
             {half_extent, -half_extent},
             {half_extent, half_extent},
             {-half_extent, half_extent}}};
}

} // namespace

TEST_CASE("shape pool creates and reads typed values", "[geometry][shape_pool]") {
    sqp::ShapePool pool;
    sqp::PhysicsMaterial circle_material;
    circle_material.base_restitution = sqp::real(0.25);
    const sqp::u32 circle_material_index = pool.materials().add(circle_material);
    const sqp::ShapeHandle circle_handle =
        pool.add_circle({4, {2, 3}, 0, 2, 1, circle_material_index});
    const sqp::ShapeHandle polygon_handle = pool.add_polygon(square(2));

    REQUIRE(pool.is_valid(circle_handle));
    REQUIRE(pool.is_valid(polygon_handle));
    CHECK(pool.kind(circle_handle) == sqp::ShapeKind::Circle);
    CHECK(pool.kind(polygon_handle) == sqp::ShapeKind::Polygon);
    CHECK(pool.shape_count() == 2);

    sqp::Circle circle;
    REQUIRE(pool.get_circle(circle_handle, circle));
    CHECK(circle.center.x == Catch::Approx(2));
    CHECK(circle.radius == Catch::Approx(4));
    CHECK_FALSE(pool.get_circle(polygon_handle, circle));

    sqp::PolygonVertices polygon;
    REQUIRE(pool.get_polygon(polygon_handle, polygon));
    REQUIRE(polygon.size() == 4);
    CHECK(polygon[0].x == Catch::Approx(-2));

    sqp::ShapeProperties properties;
    REQUIRE(pool.get_properties(polygon_handle, properties));
    CHECK(properties.area == Catch::Approx(16));
    CHECK(properties.inertia_per_unit_mass == Catch::Approx(sqp::real(8) / sqp::real(3)));
    CHECK(properties.bound_radius > properties.inner_bound_radius);
    REQUIRE(pool.get_properties(circle_handle, properties));
    CHECK(properties.density == Catch::Approx(2));
    CHECK(properties.material_index == circle_material_index);
    sqp::PhysicsMaterial actual_material;
    REQUIRE(pool.get_material(circle_handle, actual_material));
    CHECK(actual_material.base_restitution == Catch::Approx(0.25));
}

TEST_CASE("shape metadata updates do not disturb geometry", "[geometry][shape_pool]") {
    sqp::ShapePool pool;
    const sqp::ShapeHandle shape = pool.add_circle({1, {3, 4}});
    sqp::PhysicsMaterial material;
    material.dynamic_friction = sqp::real(0.125);

    REQUIRE(pool.set_density(shape, -4));
    const sqp::u32 material_index = pool.materials().add(material);
    REQUIRE(pool.set_material_index(shape, material_index));
    sqp::ShapeProperties properties;
    REQUIRE(pool.get_properties(shape, properties));
    CHECK(properties.density == 0);
    CHECK(properties.material_index == material_index);
    sqp::PhysicsMaterial actual_material;
    REQUIRE(pool.get_material(shape, actual_material));
    CHECK(actual_material.dynamic_friction == Catch::Approx(0.125));
    CHECK(properties.center_offset == sqp::Vec2(3, 4));

    REQUIRE(pool.remove(shape));
    CHECK_FALSE(pool.get_properties(shape, properties));
    CHECK_FALSE(pool.set_density(shape, 2));
    CHECK_FALSE(pool.set_material_index(shape, material_index));
}

TEST_CASE("shape pool updates transactionally and preserves type", "[geometry][shape_pool]") {
    sqp::ShapePool pool;
    const sqp::ShapeHandle circle_handle = pool.add_circle({3, {1, 2}});
    const sqp::ShapeHandle polygon_handle = pool.add_polygon(square(1));

    REQUIRE(pool.update_circle(circle_handle, {6, {4, 5}}));
    CHECK_FALSE(pool.update_circle(circle_handle, {-1, {9, 9}}));
    CHECK_FALSE(pool.update_polygon(circle_handle, square(3)));
    CHECK_FALSE(pool.update_circle(polygon_handle, {1, {}}));

    sqp::Circle circle;
    REQUIRE(pool.get_circle(circle_handle, circle));
    CHECK(circle.center == sqp::Vec2(4, 5));
    CHECK(circle.radius == Catch::Approx(6));

    const sqp::PolygonShapeDef invalid{{{0, 0}, {1, 0}}};
    CHECK_FALSE(pool.update_polygon(polygon_handle, invalid));
    sqp::PolygonVertices polygon;
    REQUIRE(pool.get_polygon(polygon_handle, polygon));
    CHECK(polygon.size() == 4);
    CHECK(polygon[0].x == Catch::Approx(-1));
}

TEST_CASE("shape pool keeps handles stable across dense removal", "[geometry][shape_pool]") {
    sqp::ShapePool pool;
    const sqp::ShapeHandle first = pool.add_circle({1, {1, 0}});
    const sqp::ShapeHandle moved = pool.add_circle({2, {2, 0}});

    REQUIRE(pool.remove(first));
    sqp::Circle circle;
    REQUIRE(pool.get_circle(moved, circle));
    CHECK(circle.radius == Catch::Approx(2));
    CHECK_FALSE(pool.is_valid(first));

    const sqp::ShapeHandle replacement = pool.add_circle({5, {}});
    CHECK(replacement.index == first.index);
    CHECK(replacement.generation != first.generation);
    CHECK(pool.is_valid(moved));

    pool.clear();
    CHECK(pool.shape_count() == 0);
    CHECK_FALSE(pool.is_valid(moved));
    CHECK_FALSE(pool.is_valid(replacement));
    const sqp::ShapeHandle after_clear = pool.add_circle({1, {}});
    CHECK(after_clear.generation != moved.generation);
}

TEST_CASE("shape pool rejects invalid fixed-capacity polygons", "[geometry][shape_pool]") {
    sqp::ShapePool pool;
    CHECK_FALSE(sqp::is_valid(pool.add_circle({0, {}})));
    CHECK_FALSE(sqp::is_valid(pool.add_polygon({{{0, 0}, {1, 0}}})));

    sqp::PolygonShapeDef too_many;
    for (std::size_t i = 0; i < sqp::config::MAX_POLYGON_VERTICES + 1; ++i) {
        too_many.vertices.push_back({static_cast<sqp::real>(i), 0});
    }
    CHECK(too_many.vertices.overflowed());
    CHECK_FALSE(sqp::is_valid(pool.add_polygon(too_many)));
    CHECK(pool.shape_count() == 0);
}

TEST_CASE("shape pool copies independently and preserves handles", "[geometry][shape_pool]") {
    STATIC_REQUIRE(std::is_copy_constructible_v<sqp::ShapePool>);
    STATIC_REQUIRE(std::is_move_constructible_v<sqp::ShapePool>);
    sqp::ShapePool original;
    const sqp::ShapeHandle handle = original.add_polygon(square(2));
    sqp::ShapePool copied = original;
    REQUIRE(copied.is_valid(handle));
    REQUIRE(original.remove(handle));
    CHECK(copied.is_valid(handle));
    sqp::ShapePool moved = std::move(original);
    CHECK_FALSE(moved.is_valid(handle));
    sqp::PolygonVertices polygon;
    REQUIRE(copied.get_polygon(handle, polygon));
    CHECK(polygon.size() == 4);
}

TEST_CASE("reserved shape pool survives dense mixed churn", "[geometry][shape_pool][stress]") {
    constexpr sqp::u32 count = 512;
    sqp::ShapePool pool;
    pool.reserve(count, count);
    std::vector<sqp::ShapeHandle> circles;
    std::vector<sqp::ShapeHandle> polygons;
    circles.reserve(count);
    polygons.reserve(count);

    for (sqp::u32 i = 0; i < count; ++i) {
        circles.push_back(pool.add_circle({sqp::real(1) + sqp::real(i % 7), {sqp::real(i), 0}}));
        sqp::PolygonShapeDef def = square(sqp::real(1) + sqp::real(i % 5));
        def.local_offset = {sqp::real(i), sqp::real(i % 11)};
        polygons.push_back(pool.add_polygon(def));
    }
    REQUIRE(pool.shape_count() == count * 2);

    bool all_removed = true;
    sqp::u32 first_failed_remove = count;
    for (sqp::u32 i = 0; i < count; i += 3) {
        const bool circle_removed = pool.remove(circles[i]);
        const bool polygon_removed = pool.remove(polygons[i]);
        const bool removed = circle_removed && polygon_removed;
        if (!removed && all_removed) {
            first_failed_remove = i;
        }
        all_removed = removed && all_removed;
    }
    CAPTURE(first_failed_remove);
    REQUIRE(all_removed);

    bool all_validity_matched = true;
    sqp::u32 first_failed_validity = count;
    for (sqp::u32 i = 0; i < count; ++i) {
        const bool expected = i % 3 != 0;
        const bool matched =
            pool.is_valid(circles[i]) == expected && pool.is_valid(polygons[i]) == expected;
        if (!matched && all_validity_matched) {
            first_failed_validity = i;
        }
        all_validity_matched = matched && all_validity_matched;
    }
    CAPTURE(first_failed_validity);
    REQUIRE(all_validity_matched);

    bool all_slots_reused = true;
    sqp::u32 first_failed_reuse = count;
    for (sqp::u32 i = 0; i < count; i += 3) {
        const sqp::ShapeHandle circle = pool.add_circle({2, {sqp::real(i), 1}});
        const sqp::ShapeHandle polygon = pool.add_polygon(square(2));
        const bool indices_valid = circle.index < count && polygon.index < count;
        const bool reused = indices_valid && circle.index % 3 == 0 && polygon.index % 3 == 0 &&
                            circle.generation != circles[circle.index].generation &&
                            polygon.generation != polygons[polygon.index].generation;
        if (!reused && all_slots_reused) {
            first_failed_reuse = i;
        }
        all_slots_reused = reused && all_slots_reused;
    }
    CAPTURE(first_failed_reuse);
    REQUIRE(all_slots_reused);
    CHECK(pool.shape_count() == count * 2);
}

namespace {
std::array<sqp::real, 7> scalar_snapshot(const sqp::ShapePool& pool, sqp::ShapeHandle h) {
    sqp::ShapeProperties value;
    REQUIRE(pool.get_properties(h, value));
    return {
        value.center_offset.x, value.center_offset.y,    value.area,   value.inertia_per_unit_mass,
        value.bound_radius,    value.inner_bound_radius, value.density};
}

void check_geometry_cache(const sqp::ShapePool& pool, sqp::ShapeHandle h) {
    sqp::ShapeProperties expected;
    if (h.kind == sqp::ShapeKind::Circle) {
        sqp::Circle circle;
        REQUIRE(pool.get_circle(h, circle));
        expected.center_offset = circle.center;
        expected.area = sqp::math::PI * circle.radius * circle.radius;
        expected.inertia_per_unit_mass = circle.radius * circle.radius / 2;
        expected.bound_radius = circle.radius;
        expected.inner_bound_radius = circle.radius;
    } else {
        sqp::PolygonVertices vertices;
        REQUIRE(pool.get_polygon(h, vertices));
        const sqp::ConvexPoly polygon(vertices);
        expected.center_offset = polygon.center_offset();
        sqp::real twice_area = 0, inertia_numerator = 0, radius_squared = 0;
        expected.inner_bound_radius = std::numeric_limits<sqp::real>::max();
        for (std::size_t i = 0; i < polygon.local.size(); ++i) {
            const auto a = polygon.local[i];
            const auto b = polygon.local[(i + 1) % polygon.local.size()];
            const auto edge = b - a;
            const auto cross = sqp::math::cross(a, b);
            twice_area += cross;
            inertia_numerator += cross * (sqp::math::dot(a, a) + sqp::math::dot(a, b) + sqp::math::dot(b, b));
            radius_squared = std::max(radius_squared, a.length_squared());
            expected.inner_bound_radius = std::min(expected.inner_bound_radius,
                                                    sqp::math::cross(edge, -a) / edge.length());
        }
        expected.area = twice_area / 2;
        expected.inertia_per_unit_mass = inertia_numerator / (6 * twice_area);
        expected.bound_radius = sqp::math::sqrt(radius_squared);
    }
    const auto actual = scalar_snapshot(pool, h);
    CHECK(actual[0] == Catch::Approx(expected.center_offset.x).margin(sqp::config::EPSILON));
    CHECK(actual[1] == Catch::Approx(expected.center_offset.y).margin(sqp::config::EPSILON));
    CHECK(actual[2] == Catch::Approx(expected.area));
    CHECK(actual[3] == Catch::Approx(expected.inertia_per_unit_mass));
    CHECK(actual[4] == Catch::Approx(expected.bound_radius));
    CHECK(actual[5] == Catch::Approx(expected.inner_bound_radius));
}
} // namespace

TEST_CASE("shape typed reads validate handles before writing outputs", "[geometry][shape_pool]") {
    sqp::ShapePool pool;
    const auto circle = pool.add_circle({2, {3, 4}, 0, 5});
    const auto polygon = pool.add_polygon(square(2));
    sqp::Circle read{{91, 92}, 93};
    REQUIRE(pool.get_circle(circle, read));
    CHECK(read.center == sqp::Vec2{3, 4});
    CHECK(read.radius == 2);
    read = {{91, 92}, 93};
    CHECK_FALSE(pool.get_circle(polygon, read));
    CHECK(read.center == sqp::Vec2{91, 92});
    CHECK(read.radius == 93);
    auto stale = circle;
    ++stale.generation;
    auto wrong_kind = circle;
    wrong_kind.kind = static_cast<sqp::ShapeKind>(255);
    bool invalid_handles_rejected = true;
    for (auto invalid : {sqp::ShapeHandle{}, stale, wrong_kind,
                         sqp::ShapeHandle{sqp::INVALID_SHAPE_INDEX, 1, sqp::ShapeKind::Circle},
                         sqp::ShapeHandle{circle.index, 0, sqp::ShapeKind::Circle}}) {
        sqp::ShapeProperties properties;
        properties.area = 91;
        invalid_handles_rejected &= !pool.get_properties(invalid, properties) &&
                                    properties.area == 91 &&
                                    !pool.set_density(invalid, sqp::real(5));
        sqp::PhysicsMaterial material;
        material.base_restitution = sqp::real(0.25);
        invalid_handles_rejected &= !pool.get_material(invalid, material) &&
                                    material.base_restitution == Catch::Approx(0.25);
    }
    CHECK(invalid_handles_rejected);
    REQUIRE(pool.remove(circle));
    CHECK_FALSE(pool.get_circle(circle, read));
    CHECK(read.radius == 93);
}

TEST_CASE("shape scalar cache follows successful geometry replacements only",
          "[geometry][shape_pool]") {
    sqp::ShapePool pool;
    auto polygon_def = square(sqp::real(2));
    polygon_def.local_offset = {5, 7};
    polygon_def.local_rotation = sqp::real(0.37);
    polygon_def.density = sqp::real(3);
    const auto polygon = pool.add_polygon(polygon_def);
    const auto circle = pool.add_circle({2, {4, 6}, 0, 7});
    check_geometry_cache(pool, circle);
    check_geometry_cache(pool, polygon);

    const auto before_circle = scalar_snapshot(pool, circle);
    const auto before_polygon = scalar_snapshot(pool, polygon);
    CHECK_FALSE(pool.update_circle(circle, {0, {99, 100}, 0, 101}));
    CHECK_FALSE(pool.update_polygon(polygon, {{{0, 0}, {1, 1}}}));
    CHECK(scalar_snapshot(pool, circle) == before_circle);
    CHECK(scalar_snapshot(pool, polygon) == before_polygon);

    polygon_def = {{{-3, -1}, {4, 0}, {2, 5}, {-2, 3}}};
    polygon_def.local_offset = {1, 2};
    polygon_def.local_rotation = sqp::real(0.19);
    polygon_def.density = sqp::real(9);
    REQUIRE(pool.update_polygon(polygon, polygon_def));
    REQUIRE(pool.update_circle(circle, {5, {-2, 9}, 0, 11}));
    CHECK(scalar_snapshot(pool, circle) != before_circle);
    CHECK(scalar_snapshot(pool, polygon) != before_polygon);
    check_geometry_cache(pool, circle);
    check_geometry_cache(pool, polygon);

    const auto h = GENERATE_COPY(circle, polygon);
    {
        auto expected = scalar_snapshot(pool, h);
        sqp::PhysicsMaterial material;
        material.dynamic_friction = sqp::real(0.13);
        const sqp::u32 material_index = pool.materials().add(material);
        REQUIRE(pool.set_material_index(h, material_index));
        CHECK(scalar_snapshot(pool, h) == expected);
        REQUIRE(pool.set_density(h, sqp::real(-1)));
        expected[6] = sqp::real(0);
        CHECK(scalar_snapshot(pool, h) == expected);
        sqp::PhysicsMaterial actual;
        REQUIRE(pool.get_material(h, actual));
        CHECK(actual.dynamic_friction == Catch::Approx(0.13));
        sqp::ShapeProperties legacy;
        REQUIRE(pool.get_properties(h, legacy));
        CHECK(legacy.area == expected[2]);
        CHECK(legacy.density == 0);
        CHECK(legacy.material_index == material_index);
    }
}

TEST_CASE("shape scalar columns and materials stay aligned through pool mutations",
          "[geometry][shape_pool]") {
    const auto kind = GENERATE(sqp::ShapeKind::Circle, sqp::ShapeKind::Polygon);
    {
        sqp::ShapePool pool;
        pool.reserve(4, 4);
        std::array<sqp::ShapeHandle, 4> handles;
        std::array<std::array<sqp::real, 7>, 4> expected;
        for (sqp::u32 i = 0; i < handles.size(); ++i) {
            sqp::PhysicsMaterial material;
            material.base_restitution = sqp::real(i + 1) / sqp::real(10);
            const sqp::u32 material_index = pool.materials().add(material);
            if (kind == sqp::ShapeKind::Circle) {
                handles[i] = pool.add_circle({sqp::real(i + 1),
                                              {sqp::real(i + 2), sqp::real(i + 3)},
                                              0,
                                              sqp::real(i + 4),
                                              1,
                                              material_index});
            } else {
                auto def = square(sqp::real(i + 1));
                def.local_offset = {sqp::real(i + 2), sqp::real(i + 3)};
                def.density = sqp::real(i + 4);
                def.material_index = material_index;
                handles[i] = pool.add_polygon(def);
            }
            expected[i] = scalar_snapshot(pool, handles[i]);
        }
        REQUIRE(pool.remove(handles[1]));
        REQUIRE(pool.remove(handles[2])); // Current tail after the previous swap-pop.
        pool.reserve(32, 32);
        bool retained_valid = true;
        for (auto i : {0u, 3u}) {
            retained_valid &= scalar_snapshot(pool, handles[i]) == expected[i];
            sqp::PhysicsMaterial material;
            retained_valid &=
                pool.get_material(handles[i], material) &&
                material.base_restitution == Catch::Approx(sqp::real(i + 1) / sqp::real(10));
        }
        CHECK(retained_valid);
        check_geometry_cache(pool, handles[0]);
        check_geometry_cache(pool, handles[3]);
        const auto replacement = kind == sqp::ShapeKind::Circle ? pool.add_circle({7, {8, 9}})
                                                                : pool.add_polygon(square(7));
        CHECK(replacement.index == handles[2].index);
        CHECK(replacement.generation != handles[2].generation);
        CHECK_FALSE(pool.is_valid(handles[2]));
        CHECK(scalar_snapshot(pool, handles[3]) == expected[3]);
        sqp::ShapePool moved = std::move(pool);
        CHECK(scalar_snapshot(moved, handles[3]) == expected[3]);
        moved.clear();
        CHECK_FALSE(moved.is_valid(handles[3]));
        const auto fresh = kind == sqp::ShapeKind::Circle ? moved.add_circle({2, {3, 4}})
                                                          : moved.add_polygon(square(2));
        check_geometry_cache(moved, fresh);
        CHECK(moved.shape_count() == 1);
    }
}

TEST_CASE("material priorities use independent full bytes", "[material]") {
    sqp::PhysicsMaterial material;
    material.friction_priority = 201;
    material.restitution_priority = 202;
    material.contact_response_priority = 203;
    material.contact_damping_priority = 204;
    CHECK(material.friction_priority == 201);
    CHECK(material.restitution_priority == 202);
    CHECK(material.contact_response_priority == 203);
    CHECK(material.contact_damping_priority == 204);
    material.friction_priority = 255;
    material.restitution_priority = 255;
    material.contact_response_priority = 255;
    material.contact_damping_priority = 255;
    CHECK(material.friction_priority == 255);
    CHECK(material.restitution_priority == 255);
    CHECK(material.contact_response_priority == 255);
    CHECK(material.contact_damping_priority == 255);
    material.friction_priority = 0;
    material.restitution_priority = 0;
    material.contact_response_priority = 0;
    material.contact_damping_priority = 0;
    CHECK(material.friction_priority == 0);
    CHECK(material.restitution_priority == 0);
    CHECK(material.contact_response_priority == 0);
    CHECK(material.contact_damping_priority == 0);
}

TEST_CASE("shape material indices share pool updates", "[material][shape_pool]") {
    sqp::ShapePool pool;
    CHECK(pool.materials().size() == 1);
    sqp::PhysicsMaterial material;
    REQUIRE(pool.materials().get(0, material));
    CHECK(material.base_restitution == 1);
    material.base_restitution = sqp::real(0.25);
    const sqp::u32 index = pool.materials().add(material);
    REQUIRE(index == 1);
    const sqp::ShapeHandle a = pool.add_circle({1, {}, 0, 1, 2, index});
    const sqp::ShapeHandle b = pool.add_circle({2, {}, 0, 1, 3, index});
    auto polygon_def = square(1);
    polygon_def.material_index = index;
    const sqp::ShapeHandle polygon = pool.add_polygon(polygon_def);
    REQUIRE(pool.is_valid(a));
    REQUIRE(pool.is_valid(b));
    REQUIRE(pool.is_valid(polygon));
    material.base_restitution = sqp::real(0.75);
    REQUIRE(pool.materials().set(index, material));
    sqp::PhysicsMaterial actual;
    REQUIRE(pool.get_material(a, actual));
    CHECK(actual.base_restitution == Catch::Approx(0.75));
    REQUIRE(pool.get_material(b, actual));
    CHECK(actual.base_restitution == Catch::Approx(0.75));
    REQUIRE(pool.get_material(polygon, actual));
    CHECK(actual.base_restitution == Catch::Approx(0.75));

    sqp::ShapeProperties properties;
    REQUIRE(pool.get_properties(a, properties));
    CHECK(properties.material_index == index);
    CHECK(properties.restitution_impact_scale == 2);
    REQUIRE(pool.set_restitution_impact_scale(a, 0));
    REQUIRE(pool.get_properties(a, properties));
    CHECK(properties.restitution_impact_scale ==
          Catch::Approx(sqp::real(1) - sqp::config::EPSILON));
    REQUIRE(pool.get_properties(b, properties));
    CHECK(properties.restitution_impact_scale == 3);
    REQUIRE(pool.get_properties(polygon, properties));
    CHECK(properties.restitution_impact_scale == 1);

    const sqp::u32 invalid = pool.materials().size();
    CHECK_FALSE(pool.set_material_index(a, invalid));
    CHECK_FALSE(pool.materials().set(invalid, material));
    CHECK_FALSE(pool.materials().get(invalid, actual));
    CHECK_FALSE(pool.is_valid(pool.add_circle({1, {}, 0, 1, 1, invalid})));
    CHECK_FALSE(pool.update_circle(a, {4, {}, 0, 1, 8, invalid}));
    polygon_def.material_index = invalid;
    CHECK_FALSE(pool.update_polygon(polygon, polygon_def));
    REQUIRE(pool.get_properties(a, properties));
    CHECK(properties.material_index == index);
    CHECK(properties.restitution_impact_scale ==
          Catch::Approx(sqp::real(1) - sqp::config::EPSILON));

    REQUIRE(pool.remove(a));
    REQUIRE(pool.get_material_index(b, properties.material_index));
    CHECK(properties.material_index == index);
    sqp::ShapePool moved = std::move(pool);
    REQUIRE(moved.get_material(b, actual));
    CHECK(actual.base_restitution == Catch::Approx(0.75));
    moved.clear();
    CHECK(moved.materials().size() == 1);
    CHECK_FALSE(moved.is_valid(b));
    REQUIRE(moved.materials().get(0, actual));
    CHECK(actual.base_restitution == 1);
}
