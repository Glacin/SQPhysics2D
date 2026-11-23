#include "world/world_stages.h"

#include <catch_amalgamated.hpp>
#include <limits>

TEST_CASE("Every built-in stage reports a specific failure", "[world][stage_error]") {
    sqp::WorldState state;

    SECTION("ProcessCommand") {
        state.settings.gravity.x = std::numeric_limits<sqp::real>::quiet_NaN();
        state.commands.heap.push_back({0, 1, sqp::PhysicsCommand{sqp::SetSettings{}}});
        CHECK_FALSE(sqp::run_process_command(state));
        CHECK(state.stage_error_detail.find("gravity") != std::string::npos);
    }
    SECTION("CollectTasks") {
        state.solve.processed_islands.push_back({{}, 2, 1});
        CHECK_FALSE(sqp::run_collect_tasks(state));
        CHECK(state.stage_error_detail.find("ends before") != std::string::npos);
    }
    SECTION("BroadPhase") {
        state.broad_phase.bodies.push_back({});
        CHECK_FALSE(sqp::run_broad_phase(state));
        CHECK(state.stage_error_detail.find("inconsistent sizes") != std::string::npos);
    }
    SECTION("NarrowPhase") {
        state.broad_phase.bodies.push_back({});
        CHECK_FALSE(sqp::run_narrow_phase(state));
        CHECK(state.stage_error_detail.find("inconsistent sizes") != std::string::npos);
    }
    SECTION("ProcessToi") {
        state.pairs.pairs->push_back({});
        bool selected = false;
        CHECK_FALSE(sqp::run_process_toi(state, 0, selected));
        CHECK(state.stage_error_detail.find("inconsistent sizes") != std::string::npos);
    }
    SECTION("UpdatePairs") {
        state.pairs.pairs->push_back({});
        CHECK_FALSE(sqp::run_update_pairs(state));
        CHECK(state.stage_error_detail.find("invalid layout") != std::string::npos);
    }
    SECTION("PreSolve") {
        state.pairs.pairs->push_back({});
        CHECK_FALSE(sqp::run_pre_solve(state));
        CHECK(state.stage_error_detail.find("invalid layout") != std::string::npos);
    }
    SECTION("Solve") {
        state.settings.gravity.x = std::numeric_limits<sqp::real>::quiet_NaN();
        CHECK_FALSE(sqp::run_solve(state));
        CHECK(state.stage_error_detail.find("gravity") != std::string::npos);
    }
    SECTION("PostSolve") {
        CHECK_FALSE(sqp::run_post_solve(state));
        CHECK(state.stage_error_detail.find("unavailable") != std::string::npos);
    }
    SECTION("AdvanceInterval") {
        const auto island = state.island_pool.create();
        state.interval.microstep = true;
        state.interval.microstep_time = std::numeric_limits<sqp::tick_t>::max();
        state.settings.interval.max_truncations_per_tick = 0;
        state.solve.processed_islands.push_back({island, 0, 0});
        CHECK_FALSE(sqp::run_advance_interval(state));
        CHECK(state.stage_error_detail.find("overflow") != std::string::npos);
    }
}

TEST_CASE("World keeps only the latest step error", "[world][stage_error]") {
    sqp::World world;
    world.state().settings.gravity.x = std::numeric_limits<sqp::real>::quiet_NaN();
    world.submit_command(0, sqp::SetSettings{});
    CHECK_FALSE(world.process_selected_commands());
    CHECK(world.last_error_message().find("ProcessCommand") != std::string_view::npos);
    CHECK(world.last_error_message().find("gravity") != std::string_view::npos);

    world.state().settings.gravity.x = 0;
    CHECK(world.process_selected_commands());
    CHECK(world.last_error_message().empty());

    world.state().schedule.time = std::numeric_limits<sqp::tick_t>::max();
    CHECK_FALSE(world.step(1));
    CHECK(world.last_error_message() == "World::step: simulation time overflow");
    world.state().schedule.time = 0;
    CHECK(world.step(0));
    CHECK(world.last_error_message().empty());
}
