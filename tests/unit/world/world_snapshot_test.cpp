#include "sqphysics2d/detail/core_world.h"

#include <catch_amalgamated.hpp>
#include <limits>

namespace {
sqp::CommandResult take(sqp::World& world, sqp::CommandId id) {
    sqp::CommandResult result;
    REQUIRE(world.take_command_result(id, result));
    return result;
}

struct CallbackContext {
    sqp::World* world;
    sqp::WorldSnapshot* snapshot;
    bool called = false;
};

void reject_reentrant_snapshot(sqp::PhysicsEvent&, void* data) {
    auto& context = *static_cast<CallbackContext*>(data);
    context.called = true;
    CHECK_FALSE(context.world->capture_snapshot(*context.snapshot));
    CHECK_FALSE(context.world->restore_snapshot(*context.snapshot));
}
} // namespace

TEST_CASE("World snapshot copies pending creation commands and completed results", "[world][snapshot]") {
    sqp::World world;
    sqp::WorldSnapshot empty, saved;
    CHECK_FALSE(world.restore_snapshot(empty));
    sqp::CircleShapeDef circle;
    circle.radius = 1;
    const auto shape_id = world.create_circle_shape(0, circle);
    const auto body_id = world.create_rigid(0);
    REQUIRE(world.capture_snapshot(saved));
    sqp::WorldSnapshot copy = saved;
    REQUIRE(world.step(0));
    sqp::WorldSnapshot completed;
    REQUIRE(world.capture_snapshot(completed));
    const auto original_shape = take(world, shape_id).shape;
    const auto original_body = take(world, body_id).body;
    REQUIRE(world.restore_snapshot(copy));
    CHECK(world.queued_command_count() == 2);
    REQUIRE(world.step(0));
    CHECK(take(world, shape_id).shape == original_shape);
    CHECK(take(world, body_id).body == original_body);
    REQUIRE(world.capture_snapshot(saved));
    CHECK(world.queued_command_count() == 0);
    REQUIRE(world.restore_snapshot(copy));
    CHECK(world.queued_command_count() == 2);
    REQUIRE(world.restore_snapshot(saved));
    CHECK(world.queued_command_count() == 0);
    REQUIRE(world.restore_snapshot(completed));
    CHECK(take(world, shape_id).shape == original_shape);
    CHECK(take(world, body_id).body == original_body);
}

TEST_CASE("World snapshot resumes a prepared interval and preserves independent bodies", "[world][snapshot]") {
    sqp::World original;
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    original.set_settings(0, settings);
    sqp::RigidBodyDef definition;
    definition.linear_velocity = {2, 0};
    const auto creation = original.create_rigid(0, definition);
    REQUIRE(original.step(0));
    const auto body = take(original, creation).body;
    const auto half = sqp::config::TU_PER_TICK / 2;
    REQUIRE(original.step(half));
    REQUIRE(original.state().schedule.prepared);
    sqp::WorldSnapshot snapshot;
    REQUIRE(original.capture_snapshot(snapshot));
    sqp::World restored;
    REQUIRE(restored.restore_snapshot(snapshot));
    CHECK(restored.time_tu() == half);
    REQUIRE(original.step(half));
    REQUIRE(restored.step(half));
    sqp::BodyState a, b;
    REQUIRE(original.query().bodies().get_state(body, a));
    REQUIRE(restored.query().bodies().get_state(body, b));
    CHECK(a.position == b.position);
    REQUIRE(restored.state().body_pool.set_pose(body, {}, 0));
    REQUIRE(original.query().bodies().get_state(body, a));
    CHECK(a.position.x != 0);
}

TEST_CASE("World snapshot rebinds joint ownership and callback data", "[world][snapshot]") {
    sqp::World original;
    auto& state = original.state();
    const auto a = state.body_pool.add_rigid();
    const auto b = state.body_pool.add_rigid();
    sqp::PinJointDef def;
    def.body_a = a;
    def.body_b = b;
    const auto joint = state.joint_pool.create(state.body_pool, def);
    REQUIRE(joint.valid());
    sqp::CircleShapeDef circle;
    circle.radius = 1;
    const auto shape = state.shape_pool.add_circle(circle);
    state.narrow_phase.candidate_cache.resize(1);
    state.narrow_phase.candidate_cache[0].view = sqp::detail::ShapeQueryView(state.shape_pool, shape);
    state.solve.contact_cache.push_back({});
    sqp::WorldSnapshot snapshot;
    CallbackContext context{&original, &snapshot};
    state.callbacks.pre_solve = reject_reentrant_snapshot;
    state.callbacks.pre_data = &context;
    REQUIRE(original.capture_snapshot(snapshot));
    sqp::World restored;
    REQUIRE(restored.restore_snapshot(snapshot));
    CHECK(restored.state().callbacks.pre_data == &context);
    CHECK(restored.state().narrow_phase.candidate_cache[0].view.pool ==
          &restored.state().shape_pool);
    CHECK(restored.state().narrow_phase.candidate_cache[0].view.is_valid());
    CHECK(restored.state().solve.contact_cache.size() == 1);
    CHECK(restored.state().joint_pool.is_valid(joint));
    CHECK(restored.state().joint_pool.first_attached(restored.state().body_pool, a) == joint);
    REQUIRE(restored.state().joint_pool.destroy(restored.state().body_pool, joint));
    CHECK(original.state().joint_pool.is_valid(joint));
}

TEST_CASE("World rejects capture after a failed step and can restore prior state", "[world][snapshot]") {
    sqp::World world;
    sqp::WorldSnapshot snapshot;
    REQUIRE(world.capture_snapshot(snapshot));
    world.state().schedule.time = std::numeric_limits<sqp::tick_t>::max();
    CHECK_FALSE(world.step(1));
    CHECK_FALSE(world.capture_snapshot(snapshot));
    REQUIRE(world.restore_snapshot(snapshot));
    CHECK(world.time_tu() == 0);
    CHECK(world.step(0));
}

TEST_CASE("World rejects snapshot operations from a physics callback", "[world][snapshot]") {
    sqp::World world;
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    world.set_settings(0, settings);
    sqp::RigidBodyDef fixed, moving;
    fixed.configs.register_type = sqp::BodyConfigs::RegisterType::ImmutableStatic;
    moving.position = {0, sqp::real(1.5)};
    moving.configs.event_flags = sqp::BodyConfigs::EventFlags::PreSolve;
    const auto fixed_id = world.create_rigid(0, fixed);
    const auto moving_id = world.create_rigid(0, moving);
    REQUIRE(world.step(0));
    const auto ground = take(world, fixed_id).body;
    const auto body = take(world, moving_id).body;
    sqp::CircleShapeDef circle;
    circle.radius = 1;
    world.create_circle_shape(0, circle, ground);
    world.create_circle_shape(0, circle, body);
    sqp::WorldSnapshot snapshot;
    CallbackContext context{&world, &snapshot};
    world.set_pre_solve_callback(0, reject_reentrant_snapshot, &context);
    REQUIRE(world.step(0));
    REQUIRE(world.capture_snapshot(snapshot));
    REQUIRE(world.step(sqp::config::TU_PER_TICK));
    CHECK(context.called);
}
