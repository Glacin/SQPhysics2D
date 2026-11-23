#include "world/world_stages.h"
#include "solver/integration.h"

#include <catch_amalgamated.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>

namespace {

sqp::u32 add_circle_proxy(
    sqp::WorldState& state, sqp::IslandId island, sqp::Vec2 position, sqp::Vec2 velocity = {},
    sqp::BodyConfigs::RegisterType registration = sqp::BodyConfigs::RegisterType::Dynamic) {
    sqp::RigidBodyDef definition;
    definition.position = position;
    definition.linear_velocity = velocity;
    definition.configs.register_type = registration;
    const auto body = state.body_pool.add_rigid(definition);
    sqp::CircleShapeDef circle;
    circle.radius = 1;
    const auto shape = state.shape_pool.add_circle(circle);
    REQUIRE(state.body_pool.attach_shape(body, shape, state.shape_pool));
    REQUIRE(state.island_pool.add_body(island, body, state.body_pool));
    const auto id = static_cast<sqp::u32>(state.broad_phase.bodies.size());
    state.broad_phase.bodies.push_back(body);
    state.broad_phase.shapes.push_back(shape);
    state.broad_phase.start_aabbs.push_back({{-20, -20}, {-19, -19}});
    state.broad_phase.end_aabbs.push_back({{-20, -20}, {-19, -19}});
    state.broad_phase.proxy_handles.push_back(state.broad_phase.tree.add_proxy(
        id, {{-20, -20}, {-19, -19}}, sqp::BroadPhase::PROXY_DYNAMIC));
    REQUIRE(state.broad_phase.proxy_handles.back() >= 0);
    return id;
}

bool has_pair(const std::vector<sqp::BroadPhase::Pair>& pairs, sqp::u32 a, sqp::u32 b) {
    return std::find(pairs.begin(), pairs.end(), sqp::BroadPhase::Pair{a, b}) != pairs.end();
}

} // namespace

TEST_CASE("microstep broad phase selects microstep islands", "[world][broad_phase]") {
    sqp::WorldState state;
    const auto regular = state.island_pool.create();
    const auto micro = state.island_pool.create();
    *state.interval.islands = {regular};
    *state.interval.delayed_islands = {micro};
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto regular_index = add_circle_proxy(state, regular, {0, 0});
    const auto micro_index = add_circle_proxy(state, micro, {5, 0});

    state.interval.islands.swap(state.interval.delayed_islands);
    state.interval.microstep = true;
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.interval_proxy_indices == std::vector<sqp::u32>{micro_index});

    state.interval.islands.swap(state.interval.delayed_islands);
    state.interval.microstep = false;
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.interval_proxy_indices == std::vector<sqp::u32>{regular_index});
}

TEST_CASE("TOI reprepare queries only the islands advanced by the microstep",
          "[world][broad_phase][toi]") {
    sqp::WorldState state;
    const auto changed = state.island_pool.create();
    const auto unchanged = state.island_pool.create();
    const auto other_unchanged = state.island_pool.create();
    *state.interval.islands = {changed, unchanged, other_unchanged};
    state.interval.target_end = sqp::config::TU_PER_TICK;
    const auto changed_proxy = add_circle_proxy(state, changed, {0, 0});
    const auto unchanged_proxy = add_circle_proxy(state, unchanged, {10, 0});
    const auto other_proxy = add_circle_proxy(state, other_unchanged, {11, 0});
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.interval_proxy_indices ==
          std::vector<sqp::u32>{changed_proxy, unchanged_proxy, other_proxy});
    REQUIRE(has_pair(state.broad_phase.ordinary_candidates, unchanged_proxy, other_proxy));

    sqp::BodyState moved_state;
    REQUIRE(state.body_pool.get_state(state.broad_phase.bodies[changed_proxy], moved_state));
    moved_state.position = {10, 0};
    REQUIRE(state.body_pool.set_state(state.broad_phase.bodies[changed_proxy], moved_state));

    state.broad_phase.reprepare_after_toi = true;
    state.broad_phase.toi_changed_islands = {changed};
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.interval_proxy_indices == std::vector<sqp::u32>{changed_proxy});
    CHECK(has_pair(state.broad_phase.candidates, unchanged_proxy, other_proxy));
    CHECK(has_pair(state.broad_phase.candidates, changed_proxy, unchanged_proxy));
    CHECK(has_pair(state.broad_phase.candidates, changed_proxy, other_proxy));
    CHECK(state.broad_phase.tree.proxy_matches(
        state.broad_phase.proxy_handles[unchanged_proxy], unchanged_proxy));

    state.broad_phase.toi_changed_islands.clear();
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.interval_proxy_indices.empty());
    CHECK(has_pair(state.broad_phase.candidates, unchanged_proxy, other_proxy));
    CHECK(has_pair(state.broad_phase.candidates, changed_proxy, unchanged_proxy));
}

TEST_CASE("multiple rotating shapes keep exact endpoint and swept bounds",
          "[world][broad_phase]") {
    for (bool indexed : {false, true}) {
        sqp::WorldState state;
        state.settings.broad_phase.proxy_extension_ratio = sqp::real(0.2);
        const auto island = state.island_pool.create();
        state.interval.islands->push_back(island);
        state.interval.target_end = sqp::config::TU_PER_SECOND / 4;
        sqp::RigidBodyDef definition;
        definition.position = {2, 3};
        definition.rotation = sqp::real(0.2);
        definition.linear_velocity = {4, -2};
        definition.angular_velocity = sqp::real(2);
        const auto body = state.body_pool.add_rigid(definition);
        REQUIRE(state.island_pool.add_body(island, body, state.body_pool));

        sqp::CircleShapeDef circle;
        circle.radius = sqp::real(0.4);
        circle.local_offset = {1, 0};
        sqp::PolygonShapeDef polygon;
        polygon.vertices = {{-1, -0.2f}, {1, -0.2f}, {1, 0.2f}, {-1, 0.2f}};
        const sqp::ShapeHandle shapes[] = {state.shape_pool.add_circle(circle),
                                           state.shape_pool.add_polygon(polygon)};
        for (sqp::u32 index = 0; index < 2; ++index) {
            REQUIRE(state.body_pool.attach_shape(body, shapes[index], state.shape_pool));
            const sqp::Aabb initial{{-10, -10}, {-9, -9}};
            state.broad_phase.bodies.push_back(body);
            state.broad_phase.shapes.push_back(shapes[index]);
            state.broad_phase.start_aabbs.push_back(initial);
            state.broad_phase.end_aabbs.push_back(initial);
            state.broad_phase.proxy_handles.push_back(state.broad_phase.tree.add_proxy(index,
                                                        initial, sqp::BroadPhase::PROXY_DYNAMIC));
            REQUIRE(state.broad_phase.proxy_handles.back() >= 0);
        }
        if (indexed) {
            const auto kind = static_cast<std::size_t>(body.kind());
            const auto dense = *state.body_pool.dense_index(body);
            state.broad_phase.proxies_by_dense[kind].resize(dense + 1);
            state.broad_phase.proxies_by_dense[kind][dense] = {0, 1};
            state.broad_phase.proxy_suppressed.resize(2, 0);
            state.broad_phase.filter_configurations.resize(2);
            REQUIRE(state.body_pool.get_configuration(body,
                        state.broad_phase.filter_configurations[0]));
            state.broad_phase.filter_configurations[1] =
                state.broad_phase.filter_configurations[0];
            state.broad_phase.proxy_index_ready = true;
        }

        sqp::BodyState start;
        REQUIRE(state.body_pool.get_state(body, start));
        sqp::BodyState end = start;
        sqp::predict_body_state(start.position, start.rotation, start.linear_velocity,
                                start.angular_velocity, {}, sqp::real(0),
                                state.interval.target_end, end.position, end.rotation);
        REQUIRE(sqp::run_broad_phase(state));
        CHECK(state.broad_phase.interval_proxy_indices == std::vector<sqp::u32>{0, 1});
        for (sqp::u32 index = 0; index < 2; ++index) {
            CHECK(state.broad_phase.start_aabbs[index] ==
                  state.shape_pool.world_aabb(shapes[index], start));
            CHECK(state.broad_phase.end_aabbs[index] ==
                  state.shape_pool.world_aabb(shapes[index], end));
            CHECK(state.broad_phase.tree.proxy_matches(state.broad_phase.proxy_handles[index],
                                                        index));
            const sqp::Aabb swept = state.shape_pool.swept_world_aabb(shapes[index], start, end);
            const sqp::Vec2 center = swept.center();
            const auto contains = [&](sqp::Vec2 point) {
                std::vector<sqp::u32> hits;
                state.broad_phase.tree.query_aabb({point, point}, hits);
                return std::find(hits.begin(), hits.end(), index) != hits.end();
            };
            const sqp::real infinity = std::numeric_limits<sqp::real>::infinity();
            CHECK(contains({swept.min.x, center.y}));
            CHECK(contains({swept.max.x, center.y}));
            CHECK(contains({center.x, swept.min.y}));
            CHECK(contains({center.x, swept.max.y}));
            CHECK_FALSE(contains({std::nextafter(swept.min.x, -infinity), center.y}));
            CHECK_FALSE(contains({std::nextafter(swept.max.x, infinity), center.y}));
            CHECK_FALSE(contains({center.x, std::nextafter(swept.min.y, -infinity)}));
            CHECK_FALSE(contains({center.x, std::nextafter(swept.max.y, infinity)}));
            const auto shift = [&](sqp::real fraction) {
                const sqp::Vec2 offset{swept.extent().x * fraction, 0};
                return sqp::Aabb{swept.min + offset, swept.max + offset};
            };
            bool reinserted = false;
            REQUIRE(state.broad_phase.tree.update_proxy(state.broad_phase.proxy_handles[index],
                                                        shift(sqp::real(0.15)),
                                                        sqp::BroadPhase::PROXY_DYNAMIC,
                                                        &reinserted));
            CHECK_FALSE(reinserted);
            REQUIRE(state.broad_phase.tree.update_proxy(state.broad_phase.proxy_handles[index],
                                                        shift(sqp::real(0.25)),
                                                        sqp::BroadPhase::PROXY_DYNAMIC,
                                                        &reinserted));
            CHECK(reinserted);
        }
    }
}

TEST_CASE("broad phase applies a direct settings change to existing proxies",
          "[world][broad_phase][settings]") {
    sqp::WorldState state;
    const auto island = state.island_pool.create();
    state.interval.islands->push_back(island);
    state.interval.target_end = sqp::config::TU_PER_TICK;
    const auto index = add_circle_proxy(state, island, {0, 0});
    const auto handle = state.broad_phase.proxy_handles[index];

    state.settings.broad_phase.proxy_extension_ratio = sqp::real(0.1);
    REQUIRE(sqp::run_broad_phase(state));
    bool reinserted = false;
    REQUIRE(state.broad_phase.tree.update_proxy(handle,
                                                {{sqp::real(-0.85), -1},
                                                 {sqp::real(1.15), 1}},
                                                sqp::BroadPhase::PROXY_DYNAMIC, &reinserted));
    CHECK(reinserted);

    state.settings.broad_phase.proxy_extension_ratio = sqp::real(0.5);
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.tree.proxy_matches(handle, index));
    REQUIRE(state.broad_phase.tree.update_proxy(handle,
                                                {{sqp::real(-0.6), -1},
                                                 {sqp::real(1.4), 1}},
                                                sqp::BroadPhase::PROXY_DYNAMIC, &reinserted));
    CHECK_FALSE(reinserted);
}

TEST_CASE("interval proxies refresh AABBs and append only their overlaps", "[world][broad_phase]") {
    sqp::WorldState state;
    const auto active = state.island_pool.create();
    const auto other = state.island_pool.create();
    REQUIRE(state.island_pool.set_time(active, 0));
    REQUIRE(state.island_pool.set_time(other, 0));
    *state.interval.islands = {active};
    state.interval.target_end = sqp::config::TU_PER_SECOND;

    const auto moving = add_circle_proxy(state, active, {0, 0}, {5, 0});
    const auto same_interval = add_circle_proxy(state, active, {6, 0}, {-1, 0});
    const auto outside = add_circle_proxy(state, other, {4, 0});
    const auto outside_pair = add_circle_proxy(state, other, {4, 1});
    const sqp::Aabb outside_bounds{{3, -1}, {5, 1}};
    REQUIRE(state.broad_phase.tree.update_proxy(state.broad_phase.proxy_handles[outside],
                                                outside_bounds, sqp::BroadPhase::PROXY_DYNAMIC));
    REQUIRE(state.broad_phase.tree.update_proxy(state.broad_phase.proxy_handles[outside_pair],
                                                {{3, 0}, {5, 2}}, sqp::BroadPhase::PROXY_DYNAMIC));

    const auto old_outside_start = state.broad_phase.start_aabbs[outside];
    const auto old_outside_end = state.broad_phase.end_aabbs[outside];
    state.broad_phase.candidates.push_back({77, 88});
    REQUIRE(sqp::run_broad_phase(state));
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
    CHECK(state.stats.pending_query_nanoseconds > 0);
    CHECK(state.stats.pending_update_nanoseconds > 0);
    CHECK(state.stats.pending_nodes_visited > 0);
    CHECK(state.stats.pending_proxy_reinsertions > 0);
#endif
    CHECK(state.broad_phase.start_aabbs[moving] == sqp::Aabb{{-1, -1}, {1, 1}});
    CHECK(state.broad_phase.end_aabbs[moving] == sqp::Aabb{{4, -1}, {6, 1}});
    CHECK(state.broad_phase.start_aabbs[outside] == old_outside_start);
    CHECK(state.broad_phase.end_aabbs[outside] == old_outside_end);
    CHECK_FALSE(has_pair(state.broad_phase.candidates, 77, 88));
    CHECK(has_pair(state.broad_phase.candidates, moving, same_interval));
    CHECK(has_pair(state.broad_phase.candidates, moving, outside));
    CHECK(has_pair(state.broad_phase.candidates, moving, outside_pair));
    CHECK_FALSE(has_pair(state.broad_phase.candidates, outside, outside_pair));
    CHECK(std::count(state.broad_phase.candidates.begin(), state.broad_phase.candidates.end(),
                     sqp::BroadPhase::Pair{moving, same_interval}) == 1);
    const auto first_size = state.broad_phase.candidates.size();
    const auto indices_capacity = state.broad_phase.interval_proxy_indices.capacity();
    const auto selected_capacity = state.broad_phase.selected.capacity();
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.candidates.size() == 2 * first_size);
    CHECK(state.broad_phase.interval_proxy_indices == std::vector<sqp::u32>{moving, same_interval});
    CHECK(state.broad_phase.selected.size() == state.broad_phase.bodies.size());
    CHECK(state.broad_phase.interval_proxy_indices.capacity() == indices_capacity);
    CHECK(state.broad_phase.selected.capacity() == selected_capacity);
}

TEST_CASE("static interval proxy does not pair with static proxy", "[world][broad_phase]") {
    sqp::WorldState state;
    const auto active = state.island_pool.create();
    const auto other = state.island_pool.create();
    *state.interval.islands = {active};
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto fixed = add_circle_proxy(state, active, {0, 0}, {5, 0},
                                        sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto other_fixed = add_circle_proxy(state, other, {0, 0});
    const auto dynamic = add_circle_proxy(state, other, {0, 0});
    REQUIRE(state.broad_phase.tree.update_proxy(state.broad_phase.proxy_handles[other_fixed],
                                                {{-1, -1}, {1, 1}}, sqp::BroadPhase::PROXY_STATIC));
    REQUIRE(state.broad_phase.tree.update_proxy(state.broad_phase.proxy_handles[dynamic],
                                                {{-1, -1}, {1, 1}},
                                                sqp::BroadPhase::PROXY_DYNAMIC));
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.end_aabbs[fixed] == sqp::Aabb{{-1, -1}, {1, 1}});
    CHECK_FALSE(has_pair(state.broad_phase.candidates, fixed, other_fixed));
    CHECK(has_pair(state.broad_phase.candidates, fixed, dynamic));
}

TEST_CASE("broad phase uses rotational sweep and rejects invalid intervals",
          "[world][broad_phase]") {
    sqp::WorldState state;
    const auto active = state.island_pool.create();
    *state.interval.islands = {active};
    state.interval.target_end = sqp::config::TU_PER_SECOND / 4;
    sqp::RigidBodyDef definition;
    definition.angular_velocity = sqp::real(2) * sqp::math::PI;
    const auto body = state.body_pool.add_rigid(definition);
    sqp::PolygonShapeDef shape_def;
    shape_def.vertices = {{-2, -0.1f}, {2, -0.1f}, {2, 0.1f}, {-2, 0.1f}};
    const auto shape = state.shape_pool.add_polygon(shape_def);
    REQUIRE(state.body_pool.attach_shape(body, shape, state.shape_pool));
    REQUIRE(state.island_pool.add_body(active, body, state.body_pool));
    state.broad_phase.bodies.push_back(body);
    state.broad_phase.shapes.push_back(shape);
    state.broad_phase.start_aabbs.push_back({});
    state.broad_phase.end_aabbs.push_back({});
    state.broad_phase.proxy_handles.push_back(
        state.broad_phase.tree.add_proxy(0, {{-1, -1}, {1, 1}}));
    REQUIRE(sqp::run_broad_phase(state));
    std::vector<sqp::u32> ids;
    state.broad_phase.tree.query_aabb({{2.001f, 0}, {2.002f, 0}}, ids);
    CHECK(ids == std::vector<sqp::u32>{0});

    state.interval.target_end = 0;
    REQUIRE(state.island_pool.set_time(active, 1));
    CHECK_FALSE(sqp::run_broad_phase(state));
}

TEST_CASE("empty interval cleans stale candidates without changing unrelated tree proxies",
          "[world][broad_phase]") {
    sqp::WorldState state;
    state.broad_phase.candidates.push_back({2, 3});
    const auto handle = state.broad_phase.tree.add_proxy(4, {{0, 0}, {1, 1}});
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.candidates.empty());
    CHECK(state.broad_phase.tree.proxy_matches(handle, 4));
    CHECK(state.broad_phase.selection_epoch == 0);
}

TEST_CASE("invalid interval and selected proxy data are rejected", "[world][broad_phase]") {
    sqp::WorldState state;
    const auto active = state.island_pool.create();
    *state.interval.islands = {active};
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto id = add_circle_proxy(state, active, {0, 0});

    SECTION("stale island") {
        *state.interval.islands = {{active.index, active.generation + 1}};
        CHECK_FALSE(sqp::run_broad_phase(state));
    }
    SECTION("invalid proxy handle") {
        state.broad_phase.proxy_handles[id] = -1;
        CHECK_FALSE(sqp::run_broad_phase(state));
    }
    SECTION("invalid body") {
        state.broad_phase.bodies[id] = {};
        CHECK_FALSE(sqp::run_broad_phase(state));
    }
    SECTION("handle points at another proxy") {
        const auto other = add_circle_proxy(state, active, {4, 0});
        state.broad_phase.proxy_handles[id] = state.broad_phase.proxy_handles[other];
        CHECK_FALSE(sqp::run_broad_phase(state));
    }
    SECTION("removed proxy handle") {
        REQUIRE(state.broad_phase.tree.remove_proxy(state.broad_phase.proxy_handles[id]));
        CHECK_FALSE(sqp::run_broad_phase(state));
    }
    SECTION("parallel array lengths") {
        state.broad_phase.proxy_handles.pop_back();
        CHECK_FALSE(sqp::run_broad_phase(state));
        state.broad_phase.proxy_handles.push_back(-1);
        state.broad_phase.end_aabbs.pop_back();
        CHECK_FALSE(sqp::run_broad_phase(state));
    }
    SECTION("invalid shape") {
        state.broad_phase.shapes[id] = {};
        CHECK_FALSE(sqp::run_broad_phase(state));
    }
    SECTION("later invalid proxy leaves earlier updates in place") {
        const auto later = add_circle_proxy(state, active, {4, 0});
        state.broad_phase.shapes[later] = {};
        CHECK_FALSE(sqp::run_broad_phase(state));
        CHECK(state.broad_phase.start_aabbs[id] == sqp::Aabb{{-1, -1}, {1, 1}});
        CHECK(state.broad_phase.interval_proxy_indices == std::vector<sqp::u32>{id});
        CHECK(state.broad_phase.candidates.empty());
    }
}

TEST_CASE("lite body registration allows broad phase update", "[world][broad_phase]") {
    sqp::WorldState state;
    const auto active = state.island_pool.create();
    *state.interval.islands = {active};
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    sqp::CircleShapeDef circle;
    circle.radius = 1;
    const auto shape = state.shape_pool.add_circle(circle);
    REQUIRE(state.shape_pool.is_valid(shape));
    const auto body = state.body_pool.add_lite_circle({}, shape, state.shape_pool);
    REQUIRE(state.body_pool.is_valid(body));
    REQUIRE(state.island_pool.add_body(active, body, state.body_pool));
    state.broad_phase.bodies.push_back(body);
    state.broad_phase.shapes.push_back(shape);
    state.broad_phase.start_aabbs.push_back({});
    state.broad_phase.end_aabbs.push_back({});
    state.broad_phase.proxy_handles.push_back(
        state.broad_phase.tree.add_proxy(0, {{-1, -1}, {1, 1}}));
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.start_aabbs[0] == sqp::Aabb{{-1, -1}, {1, 1}});
    CHECK(state.broad_phase.end_aabbs[0] == sqp::Aabb{{-1, -1}, {1, 1}});
    CHECK(state.broad_phase.candidates.empty());
}

TEST_CASE("selection stamps survive repeated intervals and wraparound", "[world][broad_phase]") {
    sqp::WorldState state;
    const auto first = state.island_pool.create();
    const auto second = state.island_pool.create();
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto a = add_circle_proxy(state, first, {0, 0});
    const auto b = add_circle_proxy(state, second, {0, 0});
    REQUIRE(state.broad_phase.tree.update_proxy(
        state.broad_phase.proxy_handles[a], {{-1, -1}, {1, 1}}, sqp::BroadPhase::PROXY_DYNAMIC));
    REQUIRE(state.broad_phase.tree.update_proxy(
        state.broad_phase.proxy_handles[b], {{-1, -1}, {1, 1}}, sqp::BroadPhase::PROXY_DYNAMIC));
    *state.interval.islands = {first};
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.candidates == std::vector<sqp::BroadPhase::Pair>{{a, b}});
    state.broad_phase.candidates.clear();
    *state.interval.islands = {second};
    state.broad_phase.selection_epoch = std::numeric_limits<sqp::u32>::max();
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.selection_epoch == 1);
    CHECK(state.broad_phase.interval_proxy_indices == std::vector<sqp::u32>{b});
    CHECK(state.broad_phase.candidates == std::vector<sqp::BroadPhase::Pair>{{a, b}});
}

TEST_CASE("two selected islands emit their shared pair once", "[world][broad_phase]") {
    sqp::WorldState state;
    const auto first = state.island_pool.create();
    const auto second = state.island_pool.create();
    *state.interval.islands = {first, second};
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto a = add_circle_proxy(state, first, {0, 0});
    const auto b = add_circle_proxy(state, second, {1, 0});
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.candidates == std::vector<sqp::BroadPhase::Pair>{{a, b}});
    CHECK(state.broad_phase.interval_proxy_indices == std::vector<sqp::u32>{a, b});
}

TEST_CASE("broad phase filters new and retained candidates using current body configuration",
          "[world][broad_phase]") {
    sqp::WorldState state;
    const auto active = state.island_pool.create();
    const auto other = state.island_pool.create();
    *state.interval.islands = {active};
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto a = add_circle_proxy(state, active, {0, 0});
    const auto filtered = add_circle_proxy(state, other, {0, 0});
    const auto bullet_a = add_circle_proxy(state, other, {0, 0});
    const auto bullet_b = add_circle_proxy(state, active, {0, 0});
    REQUIRE(state.broad_phase.tree.update_proxy(state.broad_phase.proxy_handles[filtered],
                                                {{-1, -1}, {1, 1}},
                                                sqp::BroadPhase::PROXY_DYNAMIC));
    REQUIRE(state.broad_phase.tree.update_proxy(state.broad_phase.proxy_handles[bullet_a],
                                                {{-1, -1}, {1, 1}},
                                                sqp::BroadPhase::PROXY_DYNAMIC));
    const auto bullet = sqp::BodyConfigs::OverrideFlags::IsLiteBullet;
    REQUIRE(state.body_pool.set_collision_filter(state.broad_phase.bodies[filtered], {2, 2}));
    REQUIRE(state.body_pool.set_override_flags(state.broad_phase.bodies[bullet_a], bullet));
    REQUIRE(state.body_pool.set_override_flags(state.broad_phase.bodies[bullet_b], bullet));
    state.broad_phase.candidates = {{a, filtered}, {bullet_a, bullet_b}, {a, bullet_a}};

    REQUIRE(sqp::run_broad_phase(state));
    CHECK_FALSE(has_pair(state.broad_phase.candidates, a, filtered));
    CHECK_FALSE(has_pair(state.broad_phase.candidates, bullet_a, bullet_b));
    REQUIRE_FALSE(state.broad_phase.candidates.empty());
    CHECK(state.broad_phase.candidates.front() == sqp::BroadPhase::Pair{a, bullet_a});
    CHECK(has_pair(state.broad_phase.candidates, a, bullet_b));

    REQUIRE(state.body_pool.set_collision_filter(state.broad_phase.bodies[filtered], {}));
    REQUIRE(state.body_pool.set_override_flags(state.broad_phase.bodies[bullet_b],
                                               sqp::BodyConfigs::OverrideFlags(0)));
    state.broad_phase.candidates.clear();
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(has_pair(state.broad_phase.candidates, a, filtered));
    CHECK(has_pair(state.broad_phase.candidates, bullet_a, bullet_b));

    REQUIRE(state.body_pool.set_collision_filter(state.broad_phase.bodies[filtered], {2, 2}));
    REQUIRE(state.body_pool.set_override_flags(state.broad_phase.bodies[bullet_b], bullet));
    state.interval.islands->clear();
    REQUIRE(sqp::run_broad_phase(state));
    CHECK_FALSE(has_pair(state.broad_phase.candidates, a, filtered));
    CHECK_FALSE(has_pair(state.broad_phase.candidates, bullet_a, bullet_b));
}

TEST_CASE("NoCollision removes and restores proxies across selected and unselected islands",
          "[world][broad_phase]") {
    sqp::WorldState state;
    const auto active = state.island_pool.create();
    const auto other = state.island_pool.create();
    *state.interval.islands = {active};
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto a = add_circle_proxy(state, active, {0, 0});
    const auto b = add_circle_proxy(state, other, {0, 0});
    const auto no_collision = sqp::BodyConfigs::OverrideFlags::NoCollision;

    REQUIRE(state.body_pool.set_override_flags(state.broad_phase.bodies[b], no_collision));
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.proxy_handles[b] == -1);
    CHECK(state.broad_phase.proxy_suppressed[b] == 1);
    CHECK(state.broad_phase.candidates.empty());
    std::vector<sqp::u32> hits;
    state.broad_phase.tree.query_aabb({{-1, -1}, {1, 1}}, hits);
    CHECK(std::find(hits.begin(), hits.end(), b) == hits.end());

    REQUIRE(state.body_pool.set_override_flags(state.broad_phase.bodies[b],
                                               sqp::BodyConfigs::OverrideFlags(0)));
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.tree.proxy_matches(state.broad_phase.proxy_handles[b], b));
    CHECK(state.broad_phase.proxy_suppressed[b] == 0);
    CHECK(has_pair(state.broad_phase.candidates, a, b));
    CHECK(state.broad_phase.interval_proxy_indices == std::vector<sqp::u32>{a});

    REQUIRE(state.body_pool.set_override_flags(state.broad_phase.bodies[a], no_collision));
    state.interval.islands->clear();
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.proxy_handles[a] == -1);
    CHECK(state.broad_phase.candidates.empty());
    hits.clear();
    state.broad_phase.tree.query_aabb({{-1, -1}, {1, 1}}, hits);
    CHECK(std::find(hits.begin(), hits.end(), a) == hits.end());
    REQUIRE(state.body_pool.set_override_flags(state.broad_phase.bodies[a],
                                               sqp::BodyConfigs::OverrideFlags(0)));
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.tree.proxy_matches(state.broad_phase.proxy_handles[a], a));
    CHECK(state.broad_phase.candidates.empty());
}

TEST_CASE("initially suppressed proxy can enter the tree when NoCollision is cleared",
          "[world][broad_phase]") {
    sqp::WorldState state;
    const auto active = state.island_pool.create();
    *state.interval.islands = {active};
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto id = add_circle_proxy(state, active, {0, 0});
    const auto no_collision = sqp::BodyConfigs::OverrideFlags::NoCollision;
    REQUIRE(state.broad_phase.tree.remove_proxy(state.broad_phase.proxy_handles[id]));
    state.broad_phase.proxy_handles[id] = -1;
    REQUIRE(state.body_pool.set_override_flags(state.broad_phase.bodies[id], no_collision));
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.proxy_suppressed[id] == 1);
    CHECK(state.broad_phase.tree.proxy_count() == 0);
    REQUIRE(state.body_pool.set_override_flags(state.broad_phase.bodies[id],
                                               sqp::BodyConfigs::OverrideFlags(0)));
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.tree.proxy_matches(state.broad_phase.proxy_handles[id], id));
}

TEST_CASE("high speed translation covers intermediate overlaps", "[world][broad_phase]") {
    sqp::WorldState state;
    const auto active = state.island_pool.create();
    const auto other = state.island_pool.create();
    *state.interval.islands = {active};
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto mover = add_circle_proxy(state, active, {0, 0}, {100, 0});
    const auto obstacle = add_circle_proxy(state, other, {50, 0}, {},
                                           sqp::BodyConfigs::RegisterType::ImmutableStatic);
    REQUIRE(state.broad_phase.tree.update_proxy(state.broad_phase.proxy_handles[obstacle],
                                                {{49, -1}, {51, 1}},
                                                sqp::BroadPhase::PROXY_STATIC));
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.end_aabbs[mover] == sqp::Aabb{{99, -1}, {101, 1}});
    CHECK(has_pair(state.broad_phase.candidates, mover, obstacle));
    std::vector<sqp::u32> hits;
    state.broad_phase.tree.query_aabb({{50, 0}, {50, 0}}, hits);
    CHECK(std::find(hits.begin(), hits.end(), mover) != hits.end());
}

TEST_CASE("many islands retain every required candidate once", "[world][broad_phase]") {
    sqp::WorldState state;
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    std::mt19937 generator(0x51A7u);
    std::uniform_int_distribution<int> coordinate(-15, 15);
    std::vector<sqp::Aabb> bounds;
    constexpr sqp::u32 count = 300;
    for (sqp::u32 i = 0; i < count; ++i) {
        const auto island = state.island_pool.create();
        const sqp::Vec2 position{static_cast<sqp::real>(coordinate(generator)),
                                 static_cast<sqp::real>(coordinate(generator))};
        const auto id = add_circle_proxy(state, island, position);
        REQUIRE(id == i);
        const sqp::Aabb box{position - sqp::Vec2{1, 1}, position + sqp::Vec2{1, 1}};
        bounds.push_back(box);
        REQUIRE(state.broad_phase.tree.update_proxy(state.broad_phase.proxy_handles[id], box,
                                                    sqp::BroadPhase::PROXY_DYNAMIC));
        if (i % 3 == 0)
            state.interval.islands->push_back(island);
    }
    REQUIRE(sqp::run_broad_phase(state));
    auto actual = state.broad_phase.candidates;
    CHECK(std::is_sorted(actual.begin(), actual.end()));
    std::sort(actual.begin(), actual.end());
    CHECK(std::adjacent_find(actual.begin(), actual.end()) == actual.end());
    std::vector<sqp::BroadPhase::Pair> expected;
    for (sqp::u32 a = 0; a < count; ++a)
        for (sqp::u32 b = a + 1; b < count; ++b)
            if ((a % 3 == 0 || b % 3 == 0) && sqp::Aabb::overlap(bounds[a], bounds[b]))
                expected.push_back({a, b});
    CHECK(actual == expected);
    const auto first_pass = state.broad_phase.candidates;
    state.broad_phase.candidates.clear();
    REQUIRE(sqp::run_broad_phase(state));
    CHECK(state.broad_phase.candidates == first_pass);
}
