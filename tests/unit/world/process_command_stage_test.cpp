#include "sqphysics2d/detail/core_world.h"
#include "sqphysics2d/world.h"
#include "world/world_stages.h"

#include <catch_amalgamated.hpp>

namespace {
void callback(sqp::PhysicsEvent&, void*) {}

sqp::CommandResult take(sqp::World& world, sqp::CommandId id) {
    sqp::CommandResult result;
    REQUIRE(world.take_command_result(id, result));
    CHECK_FALSE(world.take_command_result(id, result));
    return result;
}
} // namespace

TEST_CASE("Command heap processes only the selected time", "[world][command]") {
    sqp::World world;
    world.state().interval.target_end = 10;
    const auto expired = world.set_pre_solve_callback(9, callback);
    const auto first = world.set_pre_solve_callback(10, callback);
    const auto second = world.set_pre_solve_callback(10, nullptr);
    const auto future = world.set_after_solve_callback(11, callback);
    REQUIRE(world.queued_command_count() == 4);
    REQUIRE(world.process_selected_commands());
    CHECK(world.queued_command_count() == 1);
    CHECK(take(world, expired).status == sqp::CommandStatus::Expired);
    CHECK(take(world, first).status == sqp::CommandStatus::Applied);
    CHECK(take(world, second).status == sqp::CommandStatus::Applied);
    CHECK(world.state().callbacks.pre_solve == nullptr);
    CHECK(world.state().callbacks.post_solve == nullptr);
    world.state().interval.target_end = 11;
    REQUIRE(world.process_selected_commands());
    CHECK(take(world, future).status == sqp::CommandStatus::Applied);
    CHECK(world.state().callbacks.post_solve == callback);
}

TEST_CASE("A timed velocity command executes without advancing a frozen island",
          "[world][command][toi]") {
    sqp::World world;
    const auto created = world.create_rigid(0);
    world.state().interval.target_end = 0;
    REQUIRE(world.process_selected_commands());
    const auto body = take(world, created).body;
    const auto island = world.state().body_pool.island_of(body);
    REQUIRE(island.has_value());
    sqp::BodyState initial;
    REQUIRE(world.state().body_pool.get_state(body, initial));
    REQUIRE(sqp::freeze_island_until_next_tick(world.state(), *island, 0));

    const sqp::tick_t command_time = sqp::config::TU_PER_TICK / 4;
    const auto command = world.set_body_velocity(command_time, body, {3, 4});
    world.state().interval.target_end = command_time;
    REQUIRE(world.process_selected_commands());
    CHECK(take(world, command).status == sqp::CommandStatus::Applied);
    sqp::BodyState after;
    REQUIRE(world.state().body_pool.get_state(body, after));
    CHECK(after.position == initial.position);
    CHECK(after.linear_velocity == sqp::Vec2{3, 4});
    CHECK(world.state().island_pool.get(*island)->time == sqp::config::TU_PER_TICK);
    CHECK(world.state().island_pool.frozen_until(*island) == sqp::config::TU_PER_TICK);
}

TEST_CASE("Body and shape commands keep proxies aligned", "[world][command]") {
    sqp::World world;
    world.state().interval.target_end = 0;
    const auto body_id = world.create_rigid(0);
    REQUIRE(world.process_selected_commands());
    const auto body_result = take(world, body_id);
    REQUIRE(body_result.status == sqp::CommandStatus::Applied);
    const auto body = body_result.body;
    REQUIRE(world.query().bodies().is_valid(body));
    REQUIRE(world.query().bodies().island_of(body).has_value());

    sqp::CircleShapeDef circle;
    circle.radius = 1;
    const auto shape_id = world.create_circle_shape(0, circle, body);
    REQUIRE(world.process_selected_commands());
    const auto shape_result = take(world, shape_id);
    REQUIRE(shape_result.status == sqp::CommandStatus::Applied);
    REQUIRE(world.query().shapes().is_valid(shape_result.shape));
    REQUIRE(world.state().broad_phase.bodies.size() == 1);
    CHECK(world.state().broad_phase.bodies[0] == body);
    CHECK(world.state().broad_phase.shapes[0] == shape_result.shape);
    CHECK(world.state().broad_phase.tree.proxy_matches(world.state().broad_phase.proxy_handles[0],
                                                       0));

    const auto remove_id = world.destroy_shape(0, shape_result.shape);
    REQUIRE(world.process_selected_commands());
    CHECK(take(world, remove_id).status == sqp::CommandStatus::Applied);
    CHECK(world.state().broad_phase.bodies.empty());
    CHECK(world.state().broad_phase.tree.proxy_count() == 0);
    const auto destroy_id = world.destroy_body(0, body);
    REQUIRE(world.process_selected_commands());
    CHECK(take(world, destroy_id).status == sqp::CommandStatus::Applied);
    CHECK_FALSE(world.query().bodies().is_valid(body));
}

TEST_CASE("Lite creation requires a matching shape and impulse respects locks",
          "[world][command]") {
    sqp::World world;
    world.state().interval.microstep = true;
    world.state().interval.microstep_time = 4;
    sqp::CircleShapeDef circle;
    circle.radius = 1;
    const auto shape_id = world.create_circle_shape(4, circle);
    REQUIRE(world.process_selected_commands());
    const auto shape = take(world, shape_id).shape;
    REQUIRE(world.query().shapes().is_valid(shape));

    const auto bad = world.create_lite_circle(4, {}, {});
    sqp::LiteCircleBodyDef def;
    def.configs.override_flags = sqp::BodyConfigs::OverrideFlags::LockLinearX;
    const auto good = world.create_lite_circle(4, def, shape);
    REQUIRE(world.process_selected_commands());
    CHECK(take(world, bad).status == sqp::CommandStatus::Failed);
    const auto body = take(world, good).body;
    REQUIRE(world.query().bodies().is_valid(body));
    const auto impulse = world.apply_impulse(4, body, {5, 2});
    REQUIRE(world.process_selected_commands());
    CHECK(take(world, impulse).status == sqp::CommandStatus::Applied);
    sqp::BodyState state;
    REQUIRE(world.query().bodies().get_state(body, state));
    CHECK(state.linear_velocity.x == 0);
    CHECK(state.linear_velocity.y > 0);
}

TEST_CASE("Invalid command reports failure and later commands still run", "[world][command]") {
    sqp::World world;
    world.state().interval.target_end = 1;
    const auto invalid = world.set_body_velocity(1, {}, {1, 2});
    const auto callback_id = world.set_after_solve_callback(1, callback);
    REQUIRE(world.process_selected_commands());
    CHECK(take(world, invalid).status == sqp::CommandStatus::Failed);
    CHECK(take(world, callback_id).status == sqp::CommandStatus::Applied);
}

TEST_CASE("Command API applies body configuration force settings and callbacks",
          "[world][command]") {
    sqp::World world;
    world.state().interval.target_end = 2;
    const auto created = world.create_rigid(2);
    REQUIRE(world.process_selected_commands());
    const auto body = take(world, created).body;
    REQUIRE(world.query().bodies().is_valid(body));
    sqp::CircleShapeDef circle;
    circle.radius = 1;
    const auto shape_id = world.create_circle_shape(2, circle, body);
    REQUIRE(world.process_selected_commands());
    const auto shape = take(world, shape_id).shape;
    REQUIRE(world.query().shapes().is_valid(shape));

    const auto pose = world.set_body_pose(2, body, {3, 4}, 0.25);
    const auto velocity = world.set_body_velocity(2, body, {1, 2}, 0.5);
    const auto ccd = world.set_body_ccd_override(2, body, sqp::BodyConfigs::CcdOverride::Disabled);
    const auto flags =
        world.set_body_event_flags(2, body, sqp::BodyConfigs::EventFlags::AfterSolve);
    const auto override_flags =
        world.set_body_override_flags(2, body, sqp::BodyConfigs::OverrideFlags::NoCollision);
    const auto filter = world.set_body_collision_filter(2, body, {2, 4});
    const auto force = world.set_constant_force(2, body, {5, 6}, 1);
    const auto clear = world.clear_forces(2, body);
    const auto impulse = world.apply_impulse(2, body, {1, 0}, 0.1);
    sqp::PhysicsSettings settings;
    settings.gravity = {0, -5};
    settings.broad_phase.proxy_extension_ratio = sqp::real(0.2);
    const auto update_settings = world.set_settings(2, settings);
    const auto pre = world.set_pre_solve_callback(2, callback);
    const auto post = world.set_after_solve_callback(2, callback);
    REQUIRE(world.process_selected_commands());
    for (const auto id : {pose, velocity, ccd, flags, override_flags, filter, force, clear,
                          impulse, update_settings, pre, post})
        CHECK(take(world, id).status == sqp::CommandStatus::Applied);
    sqp::BodyState state;
    sqp::BodyConfiguration configuration;
    sqp::Vec2 stored_force;
    sqp::real torque = -1;
    REQUIRE(world.query().bodies().get_state(body, state));
    REQUIRE(world.query().bodies().get_configuration(body, configuration));
    REQUIRE(world.query().bodies().get_constant_force(body, stored_force, torque));
    CHECK(state.position == sqp::Vec2{3, 4});
    CHECK(state.linear_velocity.x > 1);
    CHECK(configuration.configs.ccd_override == sqp::BodyConfigs::CcdOverride::Disabled);
    CHECK(configuration.configs.event_flags == sqp::BodyConfigs::EventFlags::AfterSolve);
    CHECK(configuration.collision_filter.category_bits == 2);
    CHECK(stored_force == sqp::Vec2{});
    CHECK(torque == 0);
    CHECK(world.query().settings().gravity == sqp::Vec2{0, -5});
    CHECK(world.query().settings().broad_phase.proxy_extension_ratio == sqp::real(0.2));
    CHECK(world.state().callbacks.pre_solve == callback);
    CHECK(world.state().callbacks.post_solve == callback);
    CHECK(world.state().broad_phase.proxy_handles[0] == -1);
    CHECK(world.state().broad_phase.proxy_suppressed[0] == 1);
}

TEST_CASE("Proxy IDs are rebuilt after removing an earlier body", "[world][command]") {
    sqp::World world;
    world.state().interval.target_end = 7;
    const auto first_id = world.create_rigid(7);
    const auto second_id = world.create_rigid(7);
    REQUIRE(world.process_selected_commands());
    const auto first = take(world, first_id).body;
    const auto second = take(world, second_id).body;
    sqp::CircleShapeDef circle;
    circle.radius = 1;
    const auto first_shape_id = world.create_circle_shape(7, circle, first);
    const auto second_shape_id = world.create_circle_shape(7, circle, second);
    REQUIRE(world.process_selected_commands());
    const auto first_shape = take(world, first_shape_id).shape;
    const auto second_shape = take(world, second_shape_id).shape;
    REQUIRE(world.state().broad_phase.bodies.size() == 2);
    const auto removed = world.destroy_body(7, first);
    REQUIRE(world.process_selected_commands());
    CHECK(take(world, removed).status == sqp::CommandStatus::Applied);
    CHECK(world.state().broad_phase.bodies.size() == 1);
    CHECK(world.state().broad_phase.bodies[0] == second);
    CHECK(world.state().broad_phase.shapes[0] == second_shape);
    CHECK(world.state().broad_phase.tree.proxy_matches(world.state().broad_phase.proxy_handles[0],
                                                       0));
    CHECK(world.query().shapes().is_valid(first_shape));
}

TEST_CASE("Detached shapes can be attached and lite boxes use polygon geometry",
          "[world][command]") {
    sqp::World world;
    world.state().interval.target_end = 3;
    const auto rigid_id = world.create_rigid(3);
    sqp::CircleShapeDef circle;
    circle.radius = 1;
    const auto circle_id = world.create_circle_shape(3, circle);
    sqp::PolygonShapeDef polygon;
    polygon.vertices = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
    const auto polygon_id = world.create_polygon_shape(3, polygon);
    REQUIRE(world.process_selected_commands());
    const auto rigid = take(world, rigid_id).body;
    const auto detached = take(world, circle_id).shape;
    const auto box_shape = take(world, polygon_id).shape;
    REQUIRE(world.query().bodies().is_valid(rigid));
    REQUIRE(world.query().shapes().is_valid(detached));
    REQUIRE(world.query().shapes().is_valid(box_shape));
    const auto attach = world.attach_shape(3, rigid, detached);
    const auto box_id = world.create_aabb_box(3, {}, box_shape);
    REQUIRE(world.process_selected_commands());
    CHECK(take(world, attach).status == sqp::CommandStatus::Applied);
    const auto box = take(world, box_id).body;
    REQUIRE(world.query().bodies().is_valid(box));
    CHECK(world.state().broad_phase.bodies.size() == 2);
    const auto detach = world.detach_shape(3, rigid, detached);
    REQUIRE(world.process_selected_commands());
    CHECK(take(world, detach).status == sqp::CommandStatus::Applied);
    CHECK(world.state().broad_phase.bodies.size() == 1);
    CHECK(world.state().broad_phase.bodies[0] == box);
}
