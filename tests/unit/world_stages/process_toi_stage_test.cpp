#include "world/world_stages.h"

#include <catch_amalgamated.hpp>

namespace {

struct Fixture {
    sqp::WorldState state;
    sqp::IslandId island = state.island_pool.create();

    std::pair<sqp::BodyHandle, sqp::ShapeHandle> body(sqp::Vec2 position,
                                                      sqp::IslandId target = {}) {
        sqp::RigidBodyDef definition;
        definition.position = position;
        const auto handle = state.body_pool.add_rigid(definition);
        sqp::CircleShapeDef circle;
        circle.radius = 1;
        const auto shape = state.shape_pool.add_circle(circle);
        REQUIRE(state.body_pool.attach_shape(handle, shape, state.shape_pool));
        REQUIRE(
            state.island_pool.add_body(target.valid() ? target : island, handle, state.body_pool));
        return {handle, shape};
    }

    sqp::Pair pair(sqp::Vec2 a, sqp::Vec2 b, sqp::IslandId second_island = {},
                   sqp::IslandId first_island = {}) {
        const auto [body_a, shape_a] = body(a, first_island);
        const auto [body_b, shape_b] = body(b, second_island);
        return current_pair({body_a, body_b, shape_a, shape_b});
    }

    sqp::Pair current_pair(sqp::Pair pair) const {
        const auto versions = state.body_pool.motion_version_column(sqp::BodyConfigs::Kind::Rigid);
        pair.motion_version_a = versions[*state.body_pool.dense_index(pair.body_a)];
        pair.motion_version_b = versions[*state.body_pool.dense_index(pair.body_b)];
        return pair;
    }

    void push(sqp::Pair pair, sqp::tick_t time) {
        state.narrow_phase.toi_pairs.push_back(pair);
        state.narrow_phase.toi_times.push_back(time);
        state.narrow_phase.toi_normals.push_back({-1, 0});
    }
};

} // namespace

TEST_CASE("ProcessToi rejects mismatched columns before consuming records", "[world][toi]") {
    Fixture f;
    bool selected = true;
    f.state.narrow_phase.toi_pairs.push_back(f.pair({0, 0}, {1, 0}));
    CHECK_FALSE(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
    CHECK_FALSE(selected);
    CHECK(f.state.narrow_phase.toi_pairs.size() == 1);
    CHECK(f.state.process_toi.toi_heap.empty());

    f.state.narrow_phase.toi_times.push_back(0);
    f.state.narrow_phase.toi_normals.push_back({-1, 0});
    f.state.pairs.pairs->push_back(f.state.narrow_phase.toi_pairs[0]);
    CHECK_FALSE(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
    CHECK(f.state.narrow_phase.toi_pairs.size() == 1);
    CHECK(f.state.process_toi.toi_heap.empty());
}

TEST_CASE("TOI index preserves the first event and is rebuilt on copy", "[world][toi]") {
    Fixture f;
    const auto contact = f.pair({0, 0}, {1, 0});
    f.state.interval.target_end = 20;
    f.push(contact, 5);
    f.push(contact, 5);
    f.state.narrow_phase.toi_normals.back() = {0, 1};
    bool selected = false;
    REQUIRE(sqp::run_process_toi(f.state, 4, selected));
    REQUIRE(f.state.process_toi.toi_heap.size() == 1);
    CHECK(f.state.process_toi.toi_heap.front().normal == sqp::Vec2{-1, 0});
    auto copy = sqp::WorldState(f.state);
    f.push(contact, 6);
    REQUIRE(sqp::run_process_toi(f.state, 4, selected));
    REQUIRE(f.state.process_toi.toi_heap.size() == 2);

    copy.narrow_phase.toi_pairs.push_back(contact);
    copy.narrow_phase.toi_times.push_back(5);
    copy.narrow_phase.toi_normals.push_back({0, 1});
    REQUIRE(sqp::run_process_toi(copy, 4, selected));
    REQUIRE(copy.process_toi.toi_heap.size() == 1);
    CHECK(copy.process_toi.toi_heap.front().normal == sqp::Vec2{-1, 0});
    CHECK(copy.process_toi.next_sequence == 1);
}

TEST_CASE("ProcessToi enqueues records but preserves pending microstep output", "[world][toi]") {
    Fixture f;
    const auto first = f.pair({0, 0}, {1, 0});
    const auto second = f.pair({4, 0}, {5, 0});
    f.state.pairs.pairs->push_back(first);
    sqp::ContactManifold manifold;
    manifold.count = 1;
    REQUIRE(f.state.pairs.manifolds->append(manifold));
    const auto* active_pairs = f.state.pairs.pairs.get();
    const auto* active_manifolds = f.state.pairs.manifolds.get();
    f.state.interval.microstep = true;
    f.state.interval.microstep_time = 3;
    f.state.interval.target_end = 20;
    f.push(second, 5);
    bool selected = true;
    REQUIRE(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
    CHECK_FALSE(selected);
    CHECK(f.state.pairs.pairs->size() == 1);
    CHECK(f.state.pairs.manifolds->size() == 1);
    CHECK(f.state.pairs.pairs.get() == active_pairs);
    CHECK(f.state.pairs.manifolds.get() == active_manifolds);
    REQUIRE(f.state.process_toi.toi_heap.size() == 1);
    CHECK(f.state.process_toi.toi_heap.front().time == 5);
    CHECK(f.state.narrow_phase.toi_pairs.empty());
    CHECK(f.state.narrow_phase.toi_times.empty());
    CHECK(f.state.narrow_phase.toi_normals.empty());
    CHECK(f.state.interval.microstep);
    CHECK(f.state.interval.microstep_time == 3);
}

TEST_CASE("ProcessToi swaps ordinary results into delayed buffers by pointer", "[world][toi]") {
    Fixture f;
    f.state.interval.islands->push_back(f.island);
    const auto ordinary = f.pair({-4, 0}, {-3, 0});
    const auto contact = f.pair({0, 0}, {1, 0});
    f.state.pairs.pairs->push_back(ordinary);
    sqp::ContactManifold ordinary_manifold;
    ordinary_manifold.count = 1;
    ordinary_manifold.normal = {0, 1};
    REQUIRE(f.state.pairs.manifolds->append(ordinary_manifold));
    f.state.pairs.pairs->reserve(12);
    const auto capacity = f.state.pairs.pairs->capacity();
    const auto* ordinary_pairs = f.state.pairs.pairs.get();
    const auto* ordinary_manifolds = f.state.pairs.manifolds.get();
    const auto* empty_pairs = f.state.pairs.delayed_pairs.get();
    const auto* empty_manifolds = f.state.pairs.delayed_manifolds.get();
    const auto* ordinary_islands = f.state.interval.islands.get();
    const auto* empty_islands = f.state.interval.delayed_islands.get();
    f.push(contact, 3);
    f.state.interval.target_end = 3;

    bool selected = false;
    REQUIRE(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
    CHECK(selected);
    CHECK(f.state.interval.microstep);
    CHECK(f.state.pairs.pairs.get() == empty_pairs);
    CHECK(f.state.pairs.manifolds.get() == empty_manifolds);
    CHECK(f.state.pairs.delayed_pairs.get() == ordinary_pairs);
    CHECK(f.state.pairs.delayed_manifolds.get() == ordinary_manifolds);
    CHECK(f.state.interval.islands.get() == empty_islands);
    CHECK(f.state.interval.delayed_islands.get() == ordinary_islands);
    CHECK(f.state.interval.islands->empty());
    CHECK(*f.state.interval.delayed_islands == std::vector<sqp::IslandId>{f.island});
    REQUIRE(f.state.pairs.delayed_pairs->size() == 1);
    CHECK((*f.state.pairs.delayed_pairs)[0].body_a == ordinary.body_a);
    CHECK(f.state.pairs.delayed_manifolds->at(0).normal == sqp::Vec2{0, 1});
    CHECK(f.state.pairs.delayed_pairs->capacity() == capacity);
    REQUIRE(f.state.pairs.pairs->size() == 1);
    CHECK((*f.state.pairs.pairs)[0].body_a == contact.body_a);
    CHECK(f.state.pairs.manifolds->size() == 1);

    f.push(f.pair({4, 0}, {5, 0}), 4);
    f.state.interval.target_end = 4;
    REQUIRE(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
    CHECK_FALSE(selected);
    CHECK(f.state.pairs.pairs.get() == empty_pairs);
    CHECK(f.state.pairs.delayed_pairs.get() == ordinary_pairs);
    CHECK(f.state.interval.islands.get() == empty_islands);
    CHECK(f.state.interval.delayed_islands.get() == ordinary_islands);
    CHECK(f.state.pairs.pairs->size() == 1);
    CHECK(f.state.pairs.delayed_pairs->size() == 1);
}

TEST_CASE("A TOI batch only charges islands with solved impulses once", "[world][toi]") {
    Fixture f;
    f.state.settings.gravity = {};
    f.state.settings.interval.max_truncations_per_tick = 1;
    f.state.interval.target_end = sqp::config::TU_PER_TICK;
    const auto other = f.state.island_pool.create();
    const auto first = f.pair({0, 0}, {1, 0});
    const auto second = f.pair({4, 0}, {5, 0});
    const auto separating = f.pair({8, 0}, {9, 0}, other, other);
    REQUIRE(f.state.body_pool.set_velocity(first.body_a, {10, 0}, 0));
    REQUIRE(f.state.body_pool.set_velocity(second.body_a, {10, 0}, 0));
    REQUIRE(f.state.body_pool.set_velocity(separating.body_b, {10, 0}, 0));
    f.state.interval.tick = true;
    for (auto pair : {first, second, separating})
        f.push(f.current_pair(pair), 3);

    bool selected = false;
    REQUIRE(sqp::run_process_toi(f.state, 3, selected));
    REQUIRE(selected);
    CHECK_FALSE(f.state.interval.tick);
    CHECK(f.state.pairs.pairs->size() == 3);
    CHECK(f.state.island_pool.toi_truncations(f.island) == 0);
    REQUIRE(sqp::run_pre_solve(f.state));
    REQUIRE(sqp::run_solve(f.state));
    CHECK(f.state.island_pool.toi_truncations(f.island) == 1);
    CHECK(f.state.island_pool.toi_truncations(other) == 0);
    REQUIRE(sqp::run_advance_interval(f.state));
    CHECK(f.state.island_pool.frozen_until(f.island) == sqp::config::TU_PER_TICK);
    CHECK(f.state.island_pool.frozen_until(other) == 0);
    CHECK(f.state.island_pool.get(other)->time == 3);
}

TEST_CASE("Four zero impulse TOIs leave room for a later effective TOI", "[world][toi]") {
    Fixture f;
    f.state.settings.gravity = {};
    f.state.interval.target_end = sqp::config::TU_PER_TICK;
    const auto contact = f.pair({0, 0}, {1, 0});
    REQUIRE(f.state.body_pool.set_velocity(contact.body_b, {10, 0}, 0));

    for (sqp::tick_t event = 1; event <= 5; ++event) {
        if (event == 5)
            REQUIRE(f.state.body_pool.set_velocity(contact.body_a, {20, 0}, 0));
        const auto impact = event * (sqp::config::TU_PER_TICK / 8);
        f.push(f.current_pair(contact), impact);
        bool selected = false;
        REQUIRE(sqp::run_process_toi(f.state, impact, selected));
        REQUIRE(selected);
        REQUIRE(sqp::run_pre_solve(f.state));
        REQUIRE(sqp::run_solve(f.state));
        CHECK(f.state.island_pool.toi_truncations(f.island) == (event == 5 ? 1 : 0));
        REQUIRE(sqp::run_advance_interval(f.state));
        CHECK(f.state.island_pool.frozen_until(f.island) == 0);
        CHECK(f.state.island_pool.get(f.island)->time == impact);
    }
}

TEST_CASE("A solved TOI at the limit freezes its island until the next tick", "[world][toi]") {
    Fixture f;
    f.state.settings.gravity = {};
    f.state.settings.interval.max_truncations_per_tick = 1;
    const auto contact = f.pair({0, 0}, {1, 0});
    REQUIRE(f.state.body_pool.set_velocity(contact.body_a, {10, 0}, 0));
    const sqp::tick_t impact = sqp::config::TU_PER_TICK / 4;
    f.push(f.current_pair(contact), impact);
    f.state.interval.target_end = sqp::config::TU_PER_TICK;
    bool selected = false;
    REQUIRE(sqp::run_process_toi(f.state, impact, selected));
    REQUIRE(selected);
    REQUIRE(sqp::run_pre_solve(f.state));
    REQUIRE(sqp::run_solve(f.state));
    sqp::BodyState at_impact;
    REQUIRE(f.state.body_pool.get_state(contact.body_a, at_impact));
    REQUIRE(sqp::run_advance_interval(f.state));
    CHECK(f.state.island_pool.get(f.island)->time == sqp::config::TU_PER_TICK);
    CHECK(f.state.island_pool.frozen_until(f.island) == sqp::config::TU_PER_TICK);
    sqp::BodyState frozen;
    REQUIRE(f.state.body_pool.get_state(contact.body_a, frozen));
    CHECK(frozen.position == at_impact.position);
    CHECK(frozen.linear_velocity == at_impact.linear_velocity);

    f.state.interval.target_end = sqp::config::TU_PER_TICK;
    f.state.interval.tick = true;
    f.state.connectivity.dirty = false;
    f.state.connectivity.full_rebuild = false;
    REQUIRE(sqp::run_advance_interval(f.state));
    CHECK(f.state.island_pool.frozen_until(f.island) == 0);
    CHECK(f.state.island_pool.toi_truncations(f.island) == 0);
    CHECK(f.state.connectivity.dirty);
    CHECK(f.state.connectivity.full_rebuild);
}

TEST_CASE("Frozen contact propagation follows only scanned TOI pairs", "[world][toi]") {
    Fixture f;
    const auto b_island = f.state.island_pool.create();
    const auto c_island = f.state.island_pool.create();
    const auto [a, shape_a] = f.body({0, 0});
    const auto [b, shape_b] = f.body({1, 0}, b_island);
    const auto [c, shape_c] = f.body({2, 0}, c_island);
    const auto pair = [&](sqp::BodyHandle first, sqp::ShapeHandle first_shape,
                          sqp::BodyHandle second, sqp::ShapeHandle second_shape) {
        const auto versions_a = f.state.body_pool.motion_version_column(first.kind());
        const auto versions_b = f.state.body_pool.motion_version_column(second.kind());
        return sqp::Pair{first, second, first_shape, second_shape,
                         versions_a[*f.state.body_pool.dense_index(first)],
                         versions_b[*f.state.body_pool.dense_index(second)]};
    };
    f.state.schedule.time = 3;
    f.state.interval.target_end = sqp::config::TU_PER_TICK / 2;
    REQUIRE(sqp::freeze_island_until_next_tick(f.state, f.island, 3));
    CHECK(f.state.island_pool.frozen_until(b_island) == 0);
    CHECK(f.state.island_pool.frozen_until(c_island) == 0);
    f.push(pair(a, shape_a, b, shape_b), 5);
    bool selected = true;
    REQUIRE(sqp::run_process_toi(f.state, 5, selected));
    CHECK_FALSE(selected);
    CHECK(f.state.island_pool.frozen_until(b_island) == sqp::config::TU_PER_TICK);
    CHECK(f.state.island_pool.frozen_until(c_island) == 0);

    f.push(pair(b, shape_b, c, shape_c), 6);
    REQUIRE(sqp::run_process_toi(f.state, 6, selected));
    CHECK_FALSE(selected);
    CHECK(f.state.island_pool.frozen_until(c_island) == sqp::config::TU_PER_TICK);
}

TEST_CASE("A freeze retry preserves valid TOIs already taken at the same time", "[world][toi]") {
    Fixture f;
    const auto other = f.state.island_pool.create();
    const auto unrelated = f.state.island_pool.create();
    const auto first = f.pair({10, 0}, {11, 0}, unrelated, unrelated);
    const auto propagating = f.pair({0, 0}, {1, 0}, other);
    const auto last = f.pair({14, 0}, {15, 0}, unrelated, unrelated);
    f.state.schedule.time = 3;
    f.state.schedule.prepared = true;
    f.state.interval.target_end = 8;
    *f.state.interval.islands = {f.island, other, unrelated};
    REQUIRE(sqp::freeze_island_until_next_tick(f.state, f.island, 3));
    f.push(first, 5);
    f.push(f.current_pair(propagating), 5);
    f.push(last, 5);

    bool selected = false;
    REQUIRE(sqp::run_process_toi(f.state, 5, selected));
    CHECK_FALSE(selected);
    CHECK(f.state.broad_phase.freeze_retry);
    CHECK(f.state.broad_phase.toi_changed_islands == std::vector<sqp::IslandId>{other});
    CHECK(*f.state.interval.islands == std::vector<sqp::IslandId>{f.island, unrelated});
    CHECK(f.state.pairs.pairs->empty());
    REQUIRE(f.state.process_toi.toi_heap.size() == 2);
    CHECK(f.state.process_toi.toi_heap.front().sequence == 0);

    REQUIRE(sqp::run_process_toi(f.state, 5, selected));
    REQUIRE(selected);
    CHECK(*f.state.pairs.pairs == std::vector<sqp::Pair>{first, last});
    CHECK(f.state.process_toi.next_sequence == 3);
}

TEST_CASE("A TOI without a manifold cannot propagate a freeze", "[world][toi]") {
    Fixture f;
    const auto other = f.state.island_pool.create();
    const auto contact = f.pair({0, 0}, {10, 0}, other);
    f.state.schedule.time = 3;
    f.state.schedule.prepared = true;
    f.state.interval.target_end = 6;
    REQUIRE(sqp::freeze_island_until_next_tick(f.state, f.island, 3));
    f.push(f.current_pair(contact), 5);
    bool selected = true;
    REQUIRE(sqp::run_process_toi(f.state, 5, selected));
    CHECK_FALSE(selected);
    CHECK(f.state.process_toi.toi_heap.empty());
    CHECK(f.state.island_pool.frozen_until(other) == 0);
    CHECK(f.state.island_pool.get(other)->time == 0);
    CHECK(f.state.schedule.prepared);
}

TEST_CASE("Zero-depth polygon support uses the offset vertices", "[world][toi]") {
    Fixture f;
    const auto other = f.state.island_pool.create();
    const auto [circle, shape_circle] = f.body({sqp::real(5.00005), 0}, other);
    sqp::PolygonShapeDef definition;
    definition.vertices = {{2, -1}, {4, -1}, {4, 1}, {2, 1}};
    const auto polygon = f.state.body_pool.add_rigid({});
    const auto shape_polygon = f.state.shape_pool.add_polygon(definition);
    REQUIRE(f.state.body_pool.attach_shape(polygon, shape_polygon, f.state.shape_pool));
    REQUIRE(f.state.island_pool.add_body(f.island, polygon, f.state.body_pool));
    f.state.interval.target_end = 6;
    f.push(f.current_pair({polygon, circle, shape_polygon, shape_circle}), 5);
    bool selected = false;
    REQUIRE(sqp::run_process_toi(f.state, 5, selected));
    REQUIRE(selected);
    const auto manifold = f.state.pairs.manifolds->at(0);
    REQUIRE(manifold.count == 1);
    CHECK(manifold.points[0].depth == 0);
    CHECK(manifold.points[0].point.x == Catch::Approx(4).margin(0.0001));
}

TEST_CASE("An ordinary narrow phase contact propagates a freeze", "[world][toi]") {
    Fixture f;
    const auto other = f.state.island_pool.create();
    const auto contact = f.pair({0, 0}, {1, 0}, other);
    sqp::ContactManifold manifold;
    manifold.count = 1;
    manifold.normal = {1, 0};
    REQUIRE(f.state.pairs.manifolds->append(manifold));
    f.state.pairs.pairs->push_back(contact);
    f.state.schedule.time = 3;
    f.state.interval.target_end = 6;
    REQUIRE(sqp::freeze_island_until_next_tick(f.state, f.island, 3));
    bool selected = true;
    REQUIRE(sqp::run_process_toi(f.state, 6, selected));
    CHECK_FALSE(selected);
    CHECK(f.state.island_pool.frozen_until(other) == sqp::config::TU_PER_TICK);
    CHECK_FALSE(f.state.schedule.prepared);
}

TEST_CASE("A static contact endpoint does not acquire an island freeze", "[world][toi]") {
    Fixture f;
    const auto [moving, moving_shape] = f.body({0, 0});
    sqp::RigidBodyDef definition;
    definition.position = {1, 0};
    definition.configs.register_type = sqp::BodyConfigs::RegisterType::ImmutableStatic;
    const auto fixed = f.state.body_pool.add_rigid(definition);
    const auto fixed_shape = f.state.shape_pool.add_circle({sqp::real(1)});
    REQUIRE(f.state.body_pool.attach_shape(fixed, fixed_shape, f.state.shape_pool));
    f.state.schedule.time = 3;
    f.state.interval.target_end = 6;
    REQUIRE(sqp::freeze_island_until_next_tick(f.state, f.island, 3));
    const auto versions = f.state.body_pool.motion_version_column(moving.kind());
    f.push({moving, fixed, moving_shape, fixed_shape,
            versions[*f.state.body_pool.dense_index(moving)],
            versions[*f.state.body_pool.dense_index(fixed)]}, 5);
    f.state.schedule.prepared = true;
    bool selected = true;
    REQUIRE(sqp::run_process_toi(f.state, 5, selected));
    CHECK_FALSE(selected);
    CHECK_FALSE(f.state.broad_phase.freeze_retry);
    CHECK(f.state.schedule.prepared);
    CHECK_FALSE(f.state.body_pool.island_of(fixed).has_value());
    CHECK(f.state.island_pool.frozen_until(f.island) == sqp::config::TU_PER_TICK);
}

TEST_CASE("A frozen sleeping island stays scheduled for its resume boundary", "[world][toi]") {
    Fixture f;
    f.body({0, 0});
    REQUIRE(sqp::reschedule_island(f.state, f.island));
    REQUIRE(sqp::park_island(f.state, f.island));
    REQUIRE(f.state.island_pool.is_sleeping(f.island));
    f.state.schedule.time = 3;
    REQUIRE(sqp::freeze_island_until_next_tick(f.state, f.island, 3));
    CHECK_FALSE(f.state.island_pool.is_sleeping(f.island));
    CHECK(f.state.island_pool.get(f.island)->time == sqp::config::TU_PER_TICK);
}

TEST_CASE("Contact graph rebuilding does not spread a freeze without a scan", "[world][toi]") {
    Fixture f;
    const auto b_island = f.state.island_pool.create();
    const auto c_island = f.state.island_pool.create();
    const auto [a, shape_a] = f.body({0, 0});
    const auto [b, shape_b] = f.body({1, 0}, b_island);
    const auto [c, shape_c] = f.body({2, 0}, c_island);
    f.state.contact_graph.set_contact({a, b, shape_a, shape_b, 0, 0}, true);
    f.state.contact_graph.set_contact({b, c, shape_b, shape_c, 0, 0}, true);
    f.state.schedule.time = 3;
    REQUIRE(sqp::freeze_island_until_next_tick(f.state, f.island, 3));
    REQUIRE(sqp::run_rebuild_islands(f.state, true));
    const auto current_a = f.state.body_pool.island_of(a);
    const auto current_b = f.state.body_pool.island_of(b);
    const auto current_c = f.state.body_pool.island_of(c);
    REQUIRE(current_a.has_value());
    REQUIRE(current_b.has_value());
    REQUIRE(current_c.has_value());
    CHECK(*current_a != *current_b);
    CHECK(f.state.island_pool.frozen_until(*current_a) == sqp::config::TU_PER_TICK);
    CHECK(f.state.island_pool.frozen_until(*current_b) == 0);
    CHECK(f.state.island_pool.frozen_until(*current_c) == 0);
}

TEST_CASE("Exhausted Island does not block other TOI contacts", "[world][toi]") {
    Fixture f;
    f.state.settings.interval.max_truncations_per_tick = 1;
    REQUIRE(f.state.island_pool.set_toi_truncations(f.island, 1));
    f.push(f.pair({0, 0}, {1, 0}), 3);
    const auto other = f.state.island_pool.create();
    const auto [a, shape_a] = f.body({4, 0}, other);
    const auto [b, shape_b] = f.body({5, 0}, other);
    const auto versions = f.state.body_pool.motion_version_column(sqp::BodyConfigs::Kind::Rigid);
    const sqp::Pair allowed{a,
                            b,
                            shape_a,
                            shape_b,
                            versions[*f.state.body_pool.dense_index(a)],
                            versions[*f.state.body_pool.dense_index(b)]};
    f.push(allowed, 3);
    bool selected = false;
    REQUIRE(sqp::run_process_toi(f.state, 3, selected));
    REQUIRE(selected);
    REQUIRE(f.state.pairs.pairs->size() == 1);
    CHECK(f.state.pairs.pairs->at(0) == allowed);
    CHECK(f.state.island_pool.toi_truncations(f.island) == 1);
    CHECK(f.state.island_pool.toi_truncations(other) == 0);
}

TEST_CASE("TOI without a manifold does not consume Island budget", "[world][toi]") {
    Fixture f;
    f.push(f.pair({0, 0}, {10, 0}), 3);
    f.state.narrow_phase.toi_normals.back() = {};
    bool selected = false;
    REQUIRE(sqp::run_process_toi(f.state, 3, selected));
    CHECK_FALSE(selected);
    CHECK_FALSE(f.state.interval.microstep);
    CHECK(f.state.pairs.pairs->empty());
    CHECK(f.state.process_toi.toi_heap.empty());
    CHECK(f.state.island_pool.toi_truncations(f.island) == 0);
}

TEST_CASE("ProcessToi refuses occupied delayed buffers before consuming TOI input",
          "[world][toi]") {
    Fixture f;
    const auto event = f.pair({0, 0}, {1, 0});
    f.push(event, 3);
    f.state.interval.target_end = 3;
    f.state.pairs.delayed_pairs->push_back(event);
    sqp::ContactManifold manifold;
    manifold.count = 1;
    REQUIRE(f.state.pairs.delayed_manifolds->append(manifold));
    const auto* active = f.state.pairs.pairs.get();
    bool selected = true;
    CHECK_FALSE(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
    CHECK_FALSE(selected);
    CHECK(f.state.pairs.pairs.get() == active);
    CHECK(f.state.narrow_phase.toi_pairs.size() == 1);
    CHECK(f.state.process_toi.toi_heap.empty());
}

TEST_CASE("ProcessToi refuses occupied delayed islands before consuming TOI input",
          "[world][toi]") {
    Fixture f;
    f.push(f.pair({0, 0}, {1, 0}), 3);
    f.state.interval.delayed_islands->push_back(f.island);
    bool selected = true;
    CHECK_FALSE(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
    CHECK_FALSE(selected);
    CHECK(f.state.narrow_phase.toi_pairs.size() == 1);
    CHECK(f.state.process_toi.toi_heap.empty());
}

TEST_CASE("ProcessToi discards a stale event and retains the valid event in its group",
          "[world][toi]") {
    Fixture f;
    f.state.interval.islands->push_back(f.island);
    const auto ordinary = f.pair({-4, 0}, {-3, 0});
    f.state.pairs.pairs->push_back(ordinary);
    sqp::ContactManifold manifold;
    manifold.count = 1;
    REQUIRE(f.state.pairs.manifolds->append(manifold));
    const auto* ordinary_pairs = f.state.pairs.pairs.get();
    const auto* ordinary_manifolds = f.state.pairs.manifolds.get();
    const auto* ordinary_islands = f.state.interval.islands.get();
    f.push(f.pair({0, 0}, {1, 0}), 3);
    const auto later_island = f.state.island_pool.create();
    f.push(f.pair({4, 0}, {5, 0}, later_island), 3);
    REQUIRE(f.state.island_pool.set_time(later_island, 3));
    f.state.interval.target_end = 3;

    bool selected = true;
    CHECK(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
    CHECK(selected);
    CHECK(f.state.interval.microstep);
    CHECK(f.state.pairs.delayed_pairs.get() == ordinary_pairs);
    CHECK(f.state.pairs.delayed_manifolds.get() == ordinary_manifolds);
    CHECK(f.state.interval.delayed_islands.get() == ordinary_islands);
    CHECK(*f.state.interval.delayed_islands == std::vector<sqp::IslandId>{f.island});
    REQUIRE(f.state.pairs.pairs->size() == 1);
    CHECK((*f.state.pairs.delayed_pairs)[0].body_a == ordinary.body_a);
    CHECK(f.state.pairs.manifolds->size() == 1);
    CHECK(f.state.pairs.delayed_pairs->size() == 1);
    CHECK(f.state.pairs.delayed_manifolds->size() == 1);
}

TEST_CASE("ProcessToi skips stale groups and selects all current events at one time",
          "[world][toi]") {
    Fixture f;
    const auto stale = f.pair({0, 0}, {1, 0});
    const auto first = f.pair({4, 0}, {5, 0});
    const auto second = f.pair({8, 0}, {9, 0});
    const auto later = f.pair({12, 0}, {13, 0});
    f.push(later, 11);
    f.push(second, 7);
    f.push(stale, 3);
    f.push(first, 7);
    auto versions = f.state.body_pool.motion_version_column(stale.body_a.kind());
    ++versions[*f.state.body_pool.dense_index(stale.body_a)];
    f.state.interval.target_end = 7;
    bool selected = false;
    REQUIRE(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
    CHECK(selected);
    CHECK(f.state.interval.microstep);
    CHECK(f.state.interval.microstep_time == 7);
    REQUIRE(f.state.pairs.pairs->size() == 2);
    CHECK((*f.state.pairs.pairs)[0].body_a == second.body_a);
    CHECK((*f.state.pairs.pairs)[1].body_a == first.body_a);
    CHECK(f.state.pairs.manifolds->size() == 2);
    CHECK(f.state.pairs.manifolds->at(0).count == 1);
    REQUIRE(f.state.process_toi.toi_heap.size() == 1);
    CHECK(f.state.process_toi.toi_heap.front().time == 11);
}

TEST_CASE("ProcessToi skips an empty TOI time and restores ordinary results", "[world][toi]") {
    Fixture f;
    f.state.interval.islands->push_back(f.island);
    const auto ordinary = f.pair({-4, 0}, {-3, 0});
    f.state.pairs.pairs->push_back(ordinary);
    sqp::ContactManifold ordinary_manifold;
    ordinary_manifold.count = 1;
    REQUIRE(f.state.pairs.manifolds->append(ordinary_manifold));
    const auto* ordinary_pairs = f.state.pairs.pairs.get();
    const auto* ordinary_manifolds = f.state.pairs.manifolds.get();
    const auto* ordinary_islands = f.state.interval.islands.get();
    f.push(f.pair({0, 0}, {10, 0}), 4);
    f.state.interval.target_end = 4;
    bool selected = false;
    REQUIRE(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
    CHECK_FALSE(selected);
    CHECK_FALSE(f.state.interval.microstep);
    CHECK(f.state.pairs.pairs.get() == ordinary_pairs);
    CHECK(f.state.pairs.manifolds.get() == ordinary_manifolds);
    CHECK(f.state.interval.islands.get() == ordinary_islands);
    REQUIRE(f.state.pairs.pairs->size() == 1);
    CHECK(f.state.pairs.pairs->front() == ordinary);
    CHECK(f.state.pairs.manifolds->at(0).count == 1);
    CHECK(*f.state.interval.islands == std::vector<sqp::IslandId>{f.island});
    CHECK(f.state.pairs.delayed_pairs->empty());
    CHECK(f.state.pairs.delayed_manifolds->empty());
    CHECK(f.state.pairs.toi_local_contacts.empty());
    CHECK(f.state.pairs.delayed_toi_local_contacts.empty());
    CHECK(f.state.interval.delayed_islands->empty());
    CHECK(f.state.process_toi.toi_heap.empty());
    CHECK(f.state.island_pool.toi_truncations(f.island) == 0);
}

TEST_CASE("ProcessToi skips an empty time and selects a later contact", "[world][toi]") {
    Fixture f;
    const auto missed = f.pair({0, 0}, {10, 0});
    const auto contact = f.pair({20, 0}, {21, 0});
    const auto later = f.pair({30, 0}, {31, 0});
    f.push(missed, 4);
    f.push(contact, 6);
    f.push(later, 8);
    bool selected = false;
    REQUIRE(sqp::run_process_toi(f.state, 6, selected));
    CHECK(selected);
    CHECK(f.state.interval.microstep);
    CHECK(f.state.interval.microstep_time == 6);
    REQUIRE(f.state.pairs.pairs->size() == 1);
    CHECK(f.state.pairs.pairs->front() == contact);
    CHECK(f.state.pairs.manifolds->size() == 1);
    CHECK(f.state.pairs.toi_local_contacts.size() == 1);
    REQUIRE(f.state.process_toi.toi_heap.size() == 1);
    CHECK(f.state.process_toi.toi_heap.front().time == 8);
}

TEST_CASE("ProcessToi keeps pair and manifold columns aligned when contacts fail", "[world][toi]") {
    Fixture f;
    const auto first = f.pair({0, 0}, {1, 0});
    const auto missed = f.pair({4, 0}, {14, 0});
    const auto last = f.pair({18, 0}, {19, 0});
    f.push(first, 2);
    f.push(missed, 2);
    f.push(last, 2);
    f.state.interval.target_end = 2;
    bool selected = false;
    REQUIRE(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
    CHECK(selected);
    REQUIRE(f.state.pairs.pairs->size() == 2);
    CHECK((*f.state.pairs.pairs)[0].body_a == first.body_a);
    CHECK((*f.state.pairs.pairs)[1].body_a == last.body_a);
    CHECK(f.state.pairs.manifolds->size() == 2);
    CHECK(f.state.interval.islands->empty());
}

TEST_CASE("ProcessToi leaves events beyond the sub-interval boundary in the heap", "[world][toi]") {
    Fixture f;
    f.push(f.pair({0, 0}, {1, 0}), 6);
    f.state.interval.target_end = 5;
    bool selected = true;
    REQUIRE(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
    CHECK_FALSE(selected);
    CHECK_FALSE(f.state.interval.microstep);
    REQUIRE(f.state.process_toi.toi_heap.size() == 1);
    CHECK(f.state.process_toi.toi_heap.front().time == 6);
    f.state.interval.target_end = 6;
    REQUIRE(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
    CHECK(selected);
    CHECK(f.state.interval.microstep_time == 6);
}

TEST_CASE("ProcessToi discards current events at or before either island time", "[world][toi]") {
    SECTION("equal to both island times") {
        Fixture f;
        f.push(f.pair({0, 0}, {1, 0}), 1);
        REQUIRE(f.state.island_pool.set_time(f.island, 1));
        f.state.interval.target_end = 1;
        bool selected = true;
        CHECK(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
        CHECK_FALSE(selected);
    }
    SECTION("earlier than the second island") {
        Fixture f;
        const auto second = f.state.island_pool.create();
        f.push(f.pair({0, 0}, {1, 0}, second), 2);
        REQUIRE(f.state.island_pool.set_time(second, 3));
        f.state.interval.target_end = 2;
        bool selected = true;
        CHECK(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
        CHECK_FALSE(selected);
    }
    SECTION("stale versions are discarded before island time validation") {
        Fixture f;
        const auto pair = f.pair({0, 0}, {1, 0});
        f.push(pair, 0);
        auto versions = f.state.body_pool.motion_version_column(pair.body_a.kind());
        ++versions[*f.state.body_pool.dense_index(pair.body_a)];
        f.state.interval.target_end = 0;
        bool selected = true;
        REQUIRE(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
        CHECK_FALSE(selected);
    }
}

TEST_CASE("ProcessToi creates a zero-depth manifold for a quantized circle contact",
          "[world][toi]") {
    Fixture f;
    // Keep the gap within the zero-depth fallback tolerance in both precisions.
    const auto pair = f.pair({0, 0}, {2 + sqp::config::EPSILON * 4, 0});
    f.push(pair, 1);
    f.state.interval.target_end = 1;
    bool selected = false;
    REQUIRE(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
    CHECK(selected);
    REQUIRE(f.state.pairs.pairs->size() == 1);
    REQUIRE(f.state.pairs.manifolds->size() == 1);
    const auto manifold = f.state.pairs.manifolds->at(0);
    CHECK(manifold.count == 1);
    CHECK(manifold.points[0].depth == 0);
    CHECK(manifold.normal.x == Catch::Approx(1));
    CHECK(manifold.normal.y == Catch::Approx(0));
    CHECK(manifold.points[0].key.feature_a.type == sqp::ContactFeatureType::Circle);
    CHECK(manifold.points[0].key.feature_b.type == sqp::ContactFeatureType::Circle);
}

TEST_CASE("ProcessToi keeps pair-local impact contacts aligned at large offsets",
          "[world][toi]") {
    for (const bool reverse : {false, true}) {
        Fixture f;
        const auto first = f.pair({1024, 0}, {1026, 0});
        const auto second = f.pair({1024, 4}, {1026, 4});
        for (auto pair : {first, second}) {
            if (reverse) {
                std::swap(pair.body_a, pair.body_b);
                std::swap(pair.shape_a, pair.shape_b);
                std::swap(pair.motion_version_a, pair.motion_version_b);
            }
            f.push(pair, 1);
            const sqp::real y = pair.body_a == first.body_a || pair.body_b == first.body_a
                                    ? sqp::real(0)
                                    : sqp::real(4);
            f.state.narrow_phase.toi_normals.back() = reverse ? sqp::Vec2{1, 0}
                                                                : sqp::Vec2{-1, 0};
            f.state.narrow_phase.toi_impacts.push_back(
                {{reverse ? sqp::real(1026) : sqp::real(1024), y},
                 {},
                 {reverse ? sqp::real(-2) : sqp::real(2), 0},
                 0, 0, true});
        }
        f.state.interval.target_end = 1;
        bool selected = false;
        REQUIRE(sqp::run_process_toi(f.state, 1, selected));
        REQUIRE(selected);
        REQUIRE(f.state.pairs.pairs->size() == 2);
        REQUIRE(f.state.pairs.toi_local_contacts.size() == 2);
        for (std::size_t i = 0; i < 2; ++i) {
            const auto& local = f.state.pairs.toi_local_contacts[i];
            CHECK(local.valid);
            CHECK(local.manifold.count == 1);
            CHECK(local.manifold.points[0].point.x == Catch::Approx(reverse ? -1 : 1));
            CHECK(f.state.pairs.manifolds->at(i).points[0].point.x ==
                  Catch::Approx(local.origin.x + local.manifold.points[0].point.x));
        }
    }
}

TEST_CASE("ProcessToi fallback handles offset polygon faces", "[world][toi]") {
    Fixture f;
    const auto box = [&](sqp::Vec2 position) {
        sqp::RigidBodyDef definition;
        definition.position = position;
        const auto body = f.state.body_pool.add_rigid(definition);
        sqp::PolygonShapeDef shape_def;
        shape_def.vertices = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
        const auto shape = f.state.shape_pool.add_polygon(shape_def);
        REQUIRE(f.state.body_pool.attach_shape(body, shape, f.state.shape_pool));
        REQUIRE(f.state.island_pool.add_body(f.island, body, f.state.body_pool));
        return std::pair{body, shape};
    };
    const auto [body_a, shape_a] = box({0, 0});
    const auto [body_b, shape_b] = box({2 + sqp::config::EPSILON * 4, sqp::real(0.25)});
    const auto versions = f.state.body_pool.motion_version_column(sqp::BodyConfigs::Kind::Rigid);
    const sqp::Pair pair{body_a,
                         body_b,
                         shape_a,
                         shape_b,
                         versions[*f.state.body_pool.dense_index(body_a)],
                         versions[*f.state.body_pool.dense_index(body_b)]};
    f.push(pair, 1);
    f.state.interval.target_end = 1;
    bool selected = false;
    REQUIRE(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
    CHECK(selected);
    REQUIRE(f.state.pairs.manifolds->size() == 1);
    const auto manifold = f.state.pairs.manifolds->at(0);
    CHECK(manifold.points[0].depth == 0);
    CHECK(manifold.normal.x == Catch::Approx(1));
    CHECK(manifold.points[0].key.feature_a.type == sqp::ContactFeatureType::Edge);
    CHECK(manifold.points[0].key.feature_b.type == sqp::ContactFeatureType::Edge);
}

TEST_CASE("ProcessToi preserves a near-touching circle against a long wall",
          "[world][toi]") {
    Fixture f;
    sqp::RigidBodyDef wall_def;
    wall_def.position = {0, sqp::real(0.5)};
    const auto wall = f.state.body_pool.add_rigid(wall_def);
    sqp::PolygonShapeDef wall_shape_def;
    wall_shape_def.vertices = {{-1000, -sqp::real(0.5)},
                               {1000, -sqp::real(0.5)},
                               {1000, sqp::real(0.5)},
                               {-1000, sqp::real(0.5)}};
    const auto wall_shape = f.state.shape_pool.add_polygon(wall_shape_def);
    REQUIRE(f.state.body_pool.attach_shape(wall, wall_shape, f.state.shape_pool));
    REQUIRE(f.state.island_pool.add_body(f.island, wall, f.state.body_pool));

    const auto [circle, circle_shape] = f.body({137, -1 - sqp::real(1e-6)});
    const auto versions = f.state.body_pool.motion_version_column(sqp::BodyConfigs::Kind::Rigid);
    f.push({wall,
            circle,
            wall_shape,
            circle_shape,
            versions[*f.state.body_pool.dense_index(wall)],
            versions[*f.state.body_pool.dense_index(circle)]},
           1);
    f.state.narrow_phase.toi_normals.back() = {sqp::real(1.6e-5), 1};
    f.state.interval.target_end = 1;

    bool selected = false;
    REQUIRE(sqp::run_process_toi(f.state, 1, selected));
    REQUIRE(selected);
    REQUIRE(f.state.pairs.manifolds->size() == 1);
    const auto manifold = f.state.pairs.manifolds->at(0);
    CHECK(manifold.count == 1);
    CHECK(manifold.points[0].key.feature_a.type == sqp::ContactFeatureType::Edge);
    CHECK(manifold.points[0].key.feature_b.type == sqp::ContactFeatureType::Circle);
}
