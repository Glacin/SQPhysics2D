#include "sqphysics2d/detail/core_world.h"
#include "../../../examples/benchmark/scenario_runner.h"

#include <catch_amalgamated.hpp>

TEST_CASE("World hashes simulation state but excludes performance counters", "[world][hash]") {
    sqp::World a, b;
    CHECK(a.state_hash() == b.state_hash());
    a.state().stats.intervals = 123;
    a.state().stats.stages[0].nanoseconds = 999;
    a.state().body_pool.reserve(8, 8, 8);
    a.state().shape_pool.reserve(8, 8);
    CHECK(a.state_hash() == b.state_hash());

    const auto body = a.state().body_pool.add_rigid();
    REQUIRE(sqp::is_valid(body));
    CHECK(a.state_hash() != b.state_hash());
    const auto other = b.state().body_pool.add_rigid();
    REQUIRE(body == other);
    CHECK(a.state_hash() == b.state_hash());
    REQUIRE(a.state().body_pool.set_velocity(body, {1, 2}, 0));
    CHECK(a.state_hash() != b.state_hash());
}

TEST_CASE("World records each full tick, including in a multi-tick step", "[world][hash]") {
    sqp::World consecutive, combined;
    consecutive.create_rigid(0);
    combined.create_rigid(0);
    REQUIRE(consecutive.step(0));
    REQUIRE(combined.step(0));
    const auto half = sqp::config::TU_PER_TICK / 2;
    REQUIRE(consecutive.step(half));
    CHECK(consecutive.last_step_tick_hashes().empty());
    REQUIRE(consecutive.step(half));
    REQUIRE(consecutive.last_step_tick_hashes().size() == 1);
    const auto first = consecutive.last_step_tick_hashes().front();
    CHECK(first.tick == 1);
    REQUIRE(consecutive.step(sqp::config::TU_PER_TICK));
    const auto second = consecutive.last_step_tick_hashes().front();
    CHECK(second.tick == 2);

    REQUIRE(combined.step(2 * sqp::config::TU_PER_TICK));
    const auto records = combined.last_step_tick_hashes();
    REQUIRE(records.size() == 2);
    CHECK(records[0].tick == 1);
    CHECK(records[1].tick == 2);
    CHECK(records[0].hash == first.hash);
    CHECK(records[1].hash == second.hash);
}

TEST_CASE("World snapshot replay reproduces tick hashes", "[world][hash][snapshot]") {
    sqp::World world;
    const auto creation = world.create_rigid(0);
    REQUIRE(world.step(0));
    sqp::CommandResult result;
    REQUIRE(world.take_command_result(creation, result));
    sqp::WorldSnapshot snapshot;
    REQUIRE(world.capture_snapshot(snapshot));
    REQUIRE(world.step(sqp::config::TU_PER_TICK));
    const auto expected = world.last_step_tick_hashes().front().hash;
    REQUIRE(world.restore_snapshot(snapshot));
    CHECK(world.last_step_tick_hashes().empty());
    REQUIRE(world.step(sqp::config::TU_PER_TICK));
    CHECK(world.last_step_tick_hashes().front().hash == expected);
}

TEST_CASE("Determinism comparison finds the first differing tick", "[world][hash]") {
    const std::vector<sqp::TickHash> first{{1, 11}, {2, 22}, {3, 33}};
    auto second = first;
    CHECK_FALSE(sqp::benchmark::first_hash_mismatch(first, second));
    second[1].hash = 99;
    REQUIRE(sqp::benchmark::first_hash_mismatch(first, second));
    CHECK(*sqp::benchmark::first_hash_mismatch(first, second) == 1);
    second = first;
    second.pop_back();
    CHECK(*sqp::benchmark::first_hash_mismatch(first, second) == 2);
}
