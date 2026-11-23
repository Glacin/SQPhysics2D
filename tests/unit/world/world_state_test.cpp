#include "sqphysics2d/detail/core_world.h"

#include <catch_amalgamated.hpp>

TEST_CASE("WorldState broad phase types work through the public header", "[world][broad_phase]") {
    sqp::WorldState state;
    const auto body = state.body_pool.add_rigid({});
    sqp::CircleShapeDef shape_def;
    shape_def.radius = sqp::real(1);
    const auto shape = state.shape_pool.add_circle(shape_def);

    const sqp::Aabb start{{0, 0}, {2, 2}};
    const sqp::Aabb end{{1, 0}, {3, 2}};
    REQUIRE(state.shape_pool.is_valid(shape));
    state.broad_phase.bodies.push_back(body);
    state.broad_phase.shapes.push_back(shape);
    state.broad_phase.start_aabbs.push_back(start);
    state.broad_phase.end_aabbs.push_back(end);
    state.broad_phase.proxy_handles.push_back(
        state.broad_phase.tree.add_proxy(0, start, sqp::BroadPhase::PROXY_DYNAMIC));
    state.broad_phase.tree.add_proxy(1, {{1, 1}, {3, 3}}, sqp::BroadPhase::PROXY_STATIC);
    state.broad_phase.tree.build_pairs(state.broad_phase.candidates);

    CHECK(state.broad_phase.bodies[0] == body);
    CHECK(state.broad_phase.shapes[0] == shape);
    CHECK(state.broad_phase.start_aabbs[0] == start);
    CHECK(state.broad_phase.end_aabbs[0] == end);
    CHECK(state.broad_phase.proxy_handles[0] >= 0);
    REQUIRE(state.broad_phase.candidates.size() == 1);
    CHECK(state.broad_phase.candidates[0] == sqp::BroadPhase::Pair{0, 1});
}

TEST_CASE("WorldState initializes both pairs buffer sets", "[world][pairs]") {
    sqp::WorldState state;
    REQUIRE(state.interval.islands);
    REQUIRE(state.interval.delayed_islands);
    REQUIRE(state.pairs.pairs);
    REQUIRE(state.pairs.manifolds);
    REQUIRE(state.pairs.delayed_pairs);
    REQUIRE(state.pairs.delayed_manifolds);
    CHECK(state.pairs.pairs->empty());
    CHECK(state.pairs.manifolds->empty());
    CHECK(state.pairs.delayed_pairs->empty());
    CHECK(state.pairs.delayed_manifolds->empty());
    CHECK(state.interval.islands->empty());
    CHECK(state.interval.delayed_islands->empty());
}
