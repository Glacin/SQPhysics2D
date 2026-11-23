#pragma once

#include "sqphysics2d/detail/core_world.h"
#include "sqphysics2d/world.h"

#include <catch_amalgamated.hpp>

namespace {

struct ReferenceWorld {
    sqp::World world;

    explicit ReferenceWorld(sqp::PhysicsSettings settings) {
        settings.allow_sleep = false;
        const auto command = world.set_settings(0, settings);
        REQUIRE(world.step(0));
        sqp::CommandResult result;
        REQUIRE(world.take_command_result(command, result));
        REQUIRE(result.status == sqp::CommandStatus::Applied);
    }

    sqp::CommandResult take(sqp::CommandId id) {
        sqp::CommandResult result;
        REQUIRE(world.take_command_result(id, result));
        REQUIRE(result.status == sqp::CommandStatus::Applied);
        return result;
    }

    sqp::BodyHandle rigid(const sqp::RigidBodyDef& definition, sqp::real radius = sqp::real(0.5)) {
        const auto id = world.create_rigid(world.time_tu(), definition);
        REQUIRE(world.step(0));
        const auto body = take(id).body;
        sqp::CircleShapeDef shape;
        shape.radius = radius;
        const auto shape_id = world.create_circle_shape(world.time_tu(), shape, body);
        REQUIRE(world.step(0));
        take(shape_id);
        return body;
    }

    sqp::BodyHandle box(const sqp::RigidBodyDef& definition, sqp::Vec2 half) {
        const auto id = world.create_rigid(world.time_tu(), definition);
        REQUIRE(world.step(0));
        const auto body = take(id).body;
        sqp::PolygonShapeDef shape;
        shape.vertices = {
            {-half.x, -half.y}, {half.x, -half.y}, {half.x, half.y}, {-half.x, half.y}};
        const auto shape_id = world.create_polygon_shape(world.time_tu(), shape, body);
        REQUIRE(world.step(0));
        take(shape_id);
        return body;
    }

    sqp::BodyHandle lite(sqp::Vec2 position, sqp::Vec2 velocity,
                         sqp::real radius = sqp::real(0.5)) {
        sqp::CircleShapeDef shape;
        shape.radius = radius;
        const auto shape_id = world.create_circle_shape(world.time_tu(), shape);
        REQUIRE(world.step(0));
        const auto shape_handle = take(shape_id).shape;
        sqp::LiteCircleBodyDef definition;
        definition.position = position;
        definition.linear_velocity = velocity;
        const auto body_id = world.create_lite_circle(world.time_tu(), definition, shape_handle);
        REQUIRE(world.step(0));
        return take(body_id).body;
    }

    void material(sqp::real restitution, sqp::real static_friction = 0,
                  sqp::real dynamic_friction = 0) {
        sqp::PhysicsMaterial value;
        value.base_restitution = restitution;
        value.static_friction = static_friction;
        value.dynamic_friction = dynamic_friction;
        value.restitution_flags = sqp::RESTITUTION_DISABLE_IMPACT_SCALING;
        REQUIRE(world.state().shape_pool.materials().set(0, value));
    }

    sqp::BodyState state(sqp::BodyHandle body) const {
        sqp::BodyState value;
        REQUIRE(world.query().bodies().get_state(body, value));
        return value;
    }

    void ticks(sqp::tick_t count) {
        REQUIRE(world.step(count * sqp::config::TU_PER_TICK));
    }
};

} // namespace