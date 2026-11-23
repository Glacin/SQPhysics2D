#include <array>
#include <bit>
#include <vector>

#include "sqphysics2d/detail/data_storage.h"
#include <catch_amalgamated.hpp>

namespace {
sqp::ShapeHandle circle(sqp::ShapePool& shapes, sqp::real radius = 1, sqp::real density = 1,
                        sqp::Vec2 offset = {}) {
    sqp::CircleShapeDef def;
    def.radius = radius;
    def.density = density;
    def.local_offset = offset;
    return shapes.add_circle(def);
}
sqp::ShapeHandle box(sqp::ShapePool& shapes, sqp::real hx = 1, sqp::real hy = 1,
                     sqp::real density = 1) {
    sqp::PolygonShapeDef def;
    def.vertices = {{-hx, -hy}, {hx, -hy}, {hx, hy}, {-hx, hy}};
    def.density = density;
    return shapes.add_polygon(def);
}
} // namespace

static_assert(sizeof(sqp::BodyHandle) == 2 * sizeof(sqp::u32));

TEST_CASE("BodyPool persistent forces follow dense body movement", "[body][force]") {
    sqp::BodyPool pool;
    sqp::ShapePool shapes;
    const auto first = pool.add_rigid();
    const auto second = pool.add_rigid();
    sqp::Vec2 force;
    sqp::real torque = -1;
    REQUIRE(pool.get_constant_force(second, force, torque));
    CHECK(force == sqp::Vec2{});
    CHECK(torque == 0);
    REQUIRE(pool.set_constant_force(second, {3, 4}, 5));
    REQUIRE(pool.remove(first));
    REQUIRE(pool.get_constant_force(second, force, torque));
    CHECK(force == sqp::Vec2{3, 4});
    CHECK(torque == 5);
    CHECK(pool.scalar_column(sqp::BodyConfigs::Kind::Rigid,
                             sqp::BodyScalarColumn::ForceX)[*pool.dense_index(second)] == 3);
    REQUIRE(pool.clear_forces(second));
    REQUIRE(pool.get_constant_force(second, force, torque));
    CHECK(force == sqp::Vec2{});
    CHECK(torque == 0);

    const auto shape = shapes.add_circle({sqp::real(1)});
    const auto lite = pool.add_lite_circle({}, shape, shapes);
    REQUIRE(pool.set_constant_force(lite, {1, 2}));
    CHECK_FALSE(pool.set_constant_force(lite, {9, 9}, 1));
    REQUIRE(pool.get_constant_force(lite, force, torque));
    CHECK(force == sqp::Vec2{1, 2});
    CHECK(torque == 0);
    CHECK(pool.scalar_column(sqp::BodyConfigs::Kind::LiteCircle, sqp::BodyScalarColumn::Torque)
              .empty());

    const auto box_shape = box(shapes);
    const auto aabb = pool.add_aabb_box({}, box_shape, shapes);
    REQUIRE(pool.set_constant_force(aabb, {2, 3}));
    CHECK_FALSE(pool.set_constant_force(aabb, {2, 3}, 1));
    CHECK(
        pool.scalar_column(sqp::BodyConfigs::Kind::AabbBox, sqp::BodyScalarColumn::Torque).empty());
}

TEST_CASE("BodyHandle packs kind and all 28 index bits", "[body][handle]") {
    const sqp::BodyHandle invalid;
    CHECK_FALSE(sqp::is_valid(invalid));
    CHECK(invalid.kind() == sqp::BodyConfigs::Kind::None);

    for (const auto kind : {sqp::BodyConfigs::Kind::Rigid, sqp::BodyConfigs::Kind::LiteCircle,
                            sqp::BodyConfigs::Kind::AabbBox}) {
        const sqp::BodyHandle handle{0x0fffffffu, 42, kind};
        REQUIRE(sqp::is_valid(handle));
        CHECK(handle.index() == 0x0fffffffu);
        CHECK(handle.kind() == kind);
        CHECK(handle.generation == 42);
        const auto words = std::bit_cast<std::array<sqp::u32, 2>>(handle);
        CHECK(words[0] == 42);
        CHECK(words[1] == (static_cast<sqp::u32>(kind) << 28 | 0x0fffffffu));
    }

    CHECK_FALSE(sqp::is_valid({0x10000000u, 1, sqp::BodyConfigs::Kind::Rigid}));
    CHECK_FALSE(sqp::is_valid({0, 0, sqp::BodyConfigs::Kind::Rigid}));
    CHECK_FALSE(sqp::is_valid({0, 1, sqp::BodyConfigs::Kind::None}));

    sqp::BodyPool pool;
    const auto rigid = pool.add_rigid();
    CHECK_FALSE(
        pool.is_valid({rigid.index(), rigid.generation, sqp::BodyConfigs::Kind::LiteCircle}));
}

TEST_CASE("BodyPool typed state and configuration preserve write semantics", "[body]") {
    sqp::BodyPool pool;
    auto body = pool.add_rigid(sqp::RigidBodyDef{{1, 2}, 3, {4, 5}, 6});
    sqp::BodyState state;
    REQUIRE(pool.get_state(body, state));
    CHECK(state.position == sqp::Vec2{1, 2});
    CHECK(state.rotation == 3);
    CHECK(state.linear_velocity == sqp::Vec2{4, 5});
    CHECK(state.angular_velocity == 6);
    CHECK(state.is_alive());

    state = {{7, 8}, 9, {10, 11}, 12, sqp::BodyConfigs::RuntimeFlags::IsAlive};
    REQUIRE(pool.set_state(body, state));
    sqp::BodyState read;
    REQUIRE(pool.get_state(body, read));
    CHECK(read == state);
    REQUIRE(pool.set_velocity(body, {2, 3}, 4));
    REQUIRE(pool.get_state(body, read));
    CHECK(read.linear_velocity == sqp::Vec2{2, 3});
    CHECK(read.angular_velocity == 4);

    REQUIRE(pool.set_override_flags(body, sqp::BodyConfigs::OverrideFlags::LockLinearX |
                                              sqp::BodyConfigs::OverrideFlags::LockRot));
    REQUIRE(pool.set_velocity(body, {20, 21}, 22));
    REQUIRE(pool.get_state(body, read));
    CHECK(read.linear_velocity == sqp::Vec2{0, 21});
    CHECK(read.angular_velocity == 0);
    REQUIRE(pool.set_override_flags(body, sqp::BodyConfigs::OverrideFlags::NoSleep));
    REQUIRE(pool.get_state(body, read));

    REQUIRE(pool.set_drag_multiplier(body, -1));
    REQUIRE(pool.set_ccd_override(body, sqp::BodyConfigs::CcdOverride::Shapecast));
    REQUIRE(pool.set_collision_filter(body, {sqp::u64(8), sqp::u64(16)}));
    sqp::BodyConfiguration configuration;
    REQUIRE(pool.get_configuration(body, configuration));
    CHECK(configuration.drag_multiplier == 0);
    CHECK(configuration.configs.ccd_override == sqp::BodyConfigs::CcdOverride::Shapecast);
    CHECK(configuration.configs.kind == sqp::BodyConfigs::Kind::Rigid);
    CHECK(configuration.configs.runtime_flags == read.flags);
    CHECK(configuration.collision_filter.category_bits == 8);
    CHECK(configuration.collision_filter.collision_mask == 16);

    const auto events =
        sqp::BodyConfigs::EventFlags::ContactEvents | sqp::BodyConfigs::EventFlags::AfterSolve;
    REQUIRE(pool.set_event_flags(body, events));
    sqp::BodyConfigs::EventFlags read_events{};
    REQUIRE(pool.get_event_flags(body, read_events));
    CHECK(read_events == events);
    REQUIRE(pool.get_configuration(body, configuration));
    CHECK(configuration.configs.event_flags == events);
    read.set_alive(false);
    REQUIRE(pool.set_state(body, read));
    REQUIRE(pool.get_configuration(body, configuration));
    CHECK(configuration.configs.ccd_override == sqp::BodyConfigs::CcdOverride::Shapecast);
    CHECK(configuration.configs.event_flags == events);
    CHECK(configuration.configs.runtime_flags == read.flags);
}

TEST_CASE("BodyPool typed reads leave outputs unchanged for stale handles", "[body]") {
    sqp::BodyPool pool;
    auto body = pool.add_rigid();
    REQUIRE(pool.remove(body));
    sqp::BodyState state{{9, 9}, 9, {9, 9}, 9, sqp::BodyConfigs::RuntimeFlags(9)};
    sqp::BodyProperties properties{9, 9, 9, 9, {9, 9}, 9, 9, 9, {9, 9}};
    sqp::BodyConfiguration configuration;
    configuration.drag_multiplier = 9;
    configuration.configs.ccd_override = sqp::BodyConfigs::CcdOverride::Disabled;
    configuration.collision_filter = {9, 9};
    sqp::BodyConfigs::RegisterType registration = sqp::BodyConfigs::RegisterType::Kinematic;
    CHECK_FALSE(pool.get_state(body, state));
    CHECK_FALSE(pool.get_properties(body, properties));
    CHECK_FALSE(pool.get_configuration(body, configuration));
    CHECK_FALSE(pool.get_register_type(body, registration));
    CHECK(state.position == sqp::Vec2{9, 9});
    CHECK(properties.mass == 9);
    CHECK(configuration.drag_multiplier == 9);
    CHECK(registration == sqp::BodyConfigs::RegisterType::Kinematic);
}

TEST_CASE("BodyPool lite properties zero inapplicable rigid fields", "[body]") {
    sqp::ShapePool shapes;
    auto circle_shape = circle(shapes, 2, 3);
    auto box_shape = box(shapes, 2, 4, 5);
    sqp::BodyPool pool;
    auto circle_body = pool.add_lite_circle({}, circle_shape, shapes);
    auto box_body = pool.add_aabb_box({}, box_shape, shapes);
    sqp::BodyConfiguration config;
    REQUIRE(pool.get_configuration(circle_body, config));
    CHECK(config.configs.kind == sqp::BodyConfigs::Kind::LiteCircle);
    CHECK(config.configs.register_type == sqp::BodyConfigs::RegisterType::Dynamic);
    CHECK(config.configs.runtime_flags == sqp::BodyConfigs::RuntimeFlags::IsAlive);
    REQUIRE(pool.get_configuration(box_body, config));
    CHECK(config.configs.kind == sqp::BodyConfigs::Kind::AabbBox);
    CHECK(config.configs.register_type == sqp::BodyConfigs::RegisterType::Dynamic);
    sqp::BodyProperties value;
    REQUIRE(pool.get_properties(circle_body, value));
    CHECK(value.mass == Catch::Approx(sqp::math::PI * 4 * 3));
    CHECK(value.density == 3);
    CHECK(value.radius == 2);
    CHECK(value.half_extents == sqp::Vec2{});
    CHECK(value.inertia == 0);
    REQUIRE(pool.get_properties(box_body, value));
    CHECK(value.mass == Catch::Approx(sqp::real(8) * sqp::real(4) * 5));
    CHECK(value.density == 5);
    CHECK(value.radius == 0);
    CHECK(value.half_extents == sqp::Vec2{2, 4});
    sqp::BodyConfigs::RegisterType registration{};
    REQUIRE(pool.get_register_type(circle_body, registration));
    CHECK(registration == sqp::BodyConfigs::RegisterType::Dynamic);
    registration = sqp::BodyConfigs::RegisterType::ImmutableStatic;
    REQUIRE(pool.get_register_type(box_body, registration));
    CHECK(registration == sqp::BodyConfigs::RegisterType::Dynamic);
    REQUIRE(pool.remove(circle_body));
    registration = sqp::BodyConfigs::RegisterType::Kinematic;
    CHECK_FALSE(pool.get_register_type(circle_body, registration));
    CHECK(registration == sqp::BodyConfigs::RegisterType::Kinematic);
    REQUIRE(pool.remove(box_body));
    CHECK_FALSE(pool.get_register_type(box_body, registration));
    CHECK(registration == sqp::BodyConfigs::RegisterType::Kinematic);
}

TEST_CASE("BodyPool handle generations survive dense removal", "[body]") {
    sqp::BodyPool pool;
    std::array<sqp::BodyHandle, 3> bodies{pool.add_rigid(), pool.add_rigid(), pool.add_rigid()};
    REQUIRE(pool.set_pose(bodies[2], {7, 8}, 9));
    REQUIRE(pool.remove(bodies[1]));
    CHECK(pool.is_valid(bodies[0]));
    CHECK_FALSE(pool.is_valid(bodies[1]));
    CHECK(pool.is_valid(bodies[2]));
    sqp::BodyState state;
    REQUIRE(pool.get_state(bodies[2], state));
    CHECK(state.position == sqp::Vec2{7, 8});
    auto reused = pool.add_rigid();
    CHECK(reused.index() == bodies[1].index());
    CHECK(reused.generation != bodies[1].generation);
}

TEST_CASE("Rigid mass properties follow attached shapes", "[body][shape]") {
    sqp::ShapePool shapes;
    auto a = circle(shapes, 1, 2, {-2, 0});
    auto b = circle(shapes, 1, 2, {2, 0});
    sqp::BodyPool pool;
    auto body = pool.add_rigid();
    REQUIRE(pool.attach_shape(body, a, shapes));
    REQUIRE(pool.attach_shape(body, b, shapes));
    sqp::BodyProperties properties;
    REQUIRE(pool.get_properties(body, properties));
    CHECK(properties.mass == Catch::Approx(sqp::real(4) * sqp::math::PI));
    CHECK(properties.local_center == sqp::Vec2{});
    CHECK(properties.inertia > 0);
    CHECK(properties.inv_mass == Catch::Approx(sqp::real(1) / properties.mass));
    CHECK(properties.bound_radius == Catch::Approx(3));
    sqp::BodyHandle owner;
    REQUIRE(pool.get_shape_owner(a, owner));
    CHECK(owner == body);
    std::vector<sqp::ShapeHandle> attached;
    REQUIRE(pool.get_attached_shapes(body, attached));
    CHECK(attached == std::vector<sqp::ShapeHandle>{a, b});
}

TEST_CASE("Rigid attachment arena preserves order across growth and dense removal",
          "[body][shape][arena]") {
    sqp::ShapePool shapes;
    std::vector<sqp::ShapeHandle> handles;
    for (sqp::u32 i = 0; i < 9; ++i)
        handles.push_back(circle(shapes, sqp::real(1) + sqp::real(i) / sqp::real(10)));

    sqp::BodyPool pool;
    const auto first = pool.add_rigid();
    const auto removed = pool.add_rigid();
    const auto moved = pool.add_rigid();
    REQUIRE(pool.attach_shape(first, handles.back(), shapes));
    bool attached_all = true;
    for (sqp::u32 i = 0; i < 8; ++i)
        attached_all &= pool.attach_shape(moved, handles[i], shapes);
    REQUIRE(attached_all);
    REQUIRE(pool.detach_shape(moved, handles[3], shapes));
    REQUIRE(pool.remove(removed));

    std::vector<sqp::ShapeHandle> attached;
    REQUIRE(pool.get_attached_shapes(moved, attached));
    auto expected = std::vector<sqp::ShapeHandle>{handles[0], handles[1], handles[2], handles[4],
                                                  handles[5], handles[6], handles[7]};
    CHECK(attached == expected);
    bool owners_valid = true;
    for (auto shape : expected) {
        sqp::BodyHandle owner;
        owners_valid &= pool.get_shape_owner(shape, owner) && owner == moved;
    }
    CHECK(owners_valid);

    REQUIRE(shapes.remove(handles[3]));
    const auto reused = circle(shapes);
    CHECK(reused.index == handles[3].index);
    CHECK(reused.generation != handles[3].generation);
    sqp::BodyHandle owner;
    CHECK_FALSE(pool.get_shape_owner(reused, owner));
    REQUIRE(pool.attach_shape(first, reused, shapes));
    REQUIRE(pool.get_shape_owner(reused, owner));
    CHECK(owner == first);
}

TEST_CASE("Lite geometry refreshes and invalid geometry detaches", "[body][shape]") {
    sqp::ShapePool shapes;
    auto shape = circle(shapes, 1, 2);
    sqp::BodyPool pool;
    auto body = pool.add_lite_circle({}, shape, shapes);
    sqp::CircleShapeDef changed;
    changed.radius = 3;
    changed.density = 4;
    REQUIRE(shapes.update_circle(shape, changed));
    REQUIRE(pool.refresh_properties(body, shapes));
    sqp::BodyProperties properties;
    REQUIRE(pool.get_properties(body, properties));
    CHECK(properties.radius == 3);
    CHECK(properties.mass == Catch::Approx(sqp::math::PI * 9 * 4));
    changed.local_offset = {1, 0};
    REQUIRE(shapes.update_circle(shape, changed));
    CHECK_FALSE(pool.refresh_properties(body, shapes));
    std::vector<sqp::ShapeHandle> attached;
    REQUIRE(pool.get_attached_shapes(body, attached));
    CHECK(attached.empty());
}

TEST_CASE("Body register modes control inverse mass", "[body]") {
    bool modes_valid = true;
    for (auto registration :
         {sqp::BodyConfigs::RegisterType::Dynamic, sqp::BodyConfigs::RegisterType::Kinematic,
          sqp::BodyConfigs::RegisterType::ImmutableStatic}) {
        sqp::ShapePool shapes;
        auto shape = circle(shapes);
        sqp::BodyPool pool;
        sqp::RigidBodyDef def;
        def.configs.register_type = registration;
        auto body = pool.add_rigid(def);
        if (!pool.attach_shape(body, shape, shapes)) {
            modes_valid = false;
            continue;
        }
        sqp::BodyConfigs::RegisterType read{};
        modes_valid &= pool.get_register_type(body, read) && read == registration;
        sqp::BodyProperties properties;
        modes_valid &=
            pool.get_properties(body, properties) &&
            (properties.inv_mass > 0) == (registration == sqp::BodyConfigs::RegisterType::Dynamic);
    }
    CHECK(modes_valid);
}
