#include "sqphysics2d/detail/core_world.h"
#include "sqphysics2d/world.h"
#include "world/world_stages.h"

#include <catch_amalgamated.hpp>

namespace {
sqp::CommandResult take(sqp::World& world, sqp::CommandId id) {
    sqp::CommandResult result;
    REQUIRE(world.take_command_result(id, result));
    return result;
}

sqp::BodyState body_state(const sqp::World& world, sqp::BodyHandle body) {
    sqp::BodyState result;
    REQUIRE(world.query().bodies().get_state(body, result));
    return result;
}

void on_after_solve(sqp::PhysicsEvent& event, void* context) {
    if (event.type == sqp::PhysicsEventTypes::PostSolve)
        ++*static_cast<int*>(context);
}
} // namespace

TEST_CASE("World step pauses before the next island deadline and resumes on the same grid",
          "[world][step]") {
    sqp::World world;
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    const auto configuration = world.set_settings(0, settings);
    sqp::RigidBodyDef definition;
    definition.linear_velocity = {1, 0};
    const auto creation = world.create_rigid(0, definition);
    REQUIRE(world.step(0));
    CHECK(take(world, configuration).status == sqp::CommandStatus::Applied);
    const auto body = take(world, creation).body;
    REQUIRE(world.query().bodies().is_valid(body));

    const sqp::tick_t half = sqp::config::TU_PER_TICK / 2;
    REQUIRE(world.step(half));
    CHECK(world.time_tu() == half);
    CHECK(body_state(world, body).position.x == 0);
    REQUIRE(world.step(half));
    CHECK(world.time_tu() == sqp::config::TU_PER_TICK);
    CHECK(body_state(world, body).position.x ==
          Catch::Approx(sqp::math::tu_to_seconds<sqp::real>(sqp::config::TU_PER_TICK)));
    REQUIRE(world.step(sqp::config::TU_PER_TICK));
    CHECK(world.time_tu() == sqp::config::TU_PER_TICK * 2);
    CHECK(body_state(world, body).position.x ==
          Catch::Approx(sqp::math::tu_to_seconds<sqp::real>(sqp::config::TU_PER_TICK * 2)));
}

TEST_CASE("A timed command synchronizes its island before changing velocity", "[world][step]") {
    sqp::World world;
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    world.set_settings(0, settings);
    sqp::RigidBodyDef definition;
    definition.linear_velocity = {2, 0};
    const auto creation = world.create_rigid(0, definition);
    REQUIRE(world.step(0));
    const auto body = take(world, creation).body;
    const sqp::tick_t half = sqp::config::TU_PER_TICK / 2;
    const auto change = world.set_body_velocity(half, body, {}, 0);
    REQUIRE(world.step(half));
    CHECK(take(world, change).status == sqp::CommandStatus::Applied);
    CHECK(body_state(world, body).position.x ==
          Catch::Approx(sqp::math::tu_to_seconds<sqp::real>(half) * 2));
    REQUIRE(world.step(sqp::config::TU_PER_TICK));
    CHECK(body_state(world, body).position.x ==
          Catch::Approx(sqp::math::tu_to_seconds<sqp::real>(half) * 2));
}

TEST_CASE("Island periods select only the bodies due at each grid time", "[world][step]") {
    sqp::World world;
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    world.set_settings(0, settings);
    sqp::RigidBodyDef fast, slow;
    fast.linear_velocity = {100, 0};
    slow.linear_velocity = {1, 0};
    const auto fast_id = world.create_rigid(0, fast);
    const auto slow_id = world.create_rigid(0, slow);
    REQUIRE(world.step(0));
    const auto fast_body = take(world, fast_id).body;
    const auto slow_body = take(world, slow_id).body;
    const sqp::tick_t short_period = settings.interval.min_period_tu;
    REQUIRE(world.step(short_period));
    CHECK(body_state(world, fast_body).position.x > 0);
    CHECK(body_state(world, slow_body).position.x == 0);
    REQUIRE(world.step(sqp::config::TU_PER_TICK - short_period));
    CHECK(body_state(world, slow_body).position.x > 0);
}

TEST_CASE("Bodies using CCD do not shorten ordinary island periods", "[world][step][ccd]") {
    using Mode = sqp::BodyConfigs::CcdOverride;
    const auto period_for = [](Mode mode, sqp::real speed, sqp::real turn,
                               bool enable_circle = true) {
        sqp::WorldState state;
        state.settings.gravity = {};
        state.settings.ccd.enable_circle = enable_circle;
        sqp::RigidBodyDef definition;
        definition.linear_velocity = {speed, 0};
        definition.angular_velocity = turn;
        definition.configs.ccd_override = mode;
        const auto body = state.body_pool.add_rigid(definition);
        sqp::CircleShapeDef circle;
        circle.radius = 1;
        const auto shape = state.shape_pool.add_circle(circle);
        REQUIRE(state.body_pool.attach_shape(body, shape, state.shape_pool));
        const auto island = state.island_pool.create();
        REQUIRE(state.island_pool.add_body(island, body, state.body_pool));
        REQUIRE(sqp::reschedule_island(state, island));
        return state.schedule.entries[island.index].period;
    };
    const auto tick = sqp::config::TU_PER_TICK;
    const auto short_period = tick / 4;
    CHECK(period_for(Mode::Shapecast, 100, 0) == tick);
    CHECK(period_for(Mode::Auto, 100, 0) == tick);
    CHECK(period_for(Mode::Auto, 0, 100) == tick);
    CHECK(period_for(Mode::Disabled, 100, 0) == short_period);
    CHECK(period_for(Mode::Auto, 10, 0) == tick / 2);
    CHECK(period_for(Mode::Shapecast, 100, 0, false) == short_period);

    sqp::WorldState mixed;
    mixed.settings.gravity = {};
    const auto island = mixed.island_pool.create();
    sqp::CircleShapeDef circle;
    circle.radius = 1;
    for (Mode mode : {Mode::Shapecast, Mode::Disabled}) {
        sqp::RigidBodyDef definition;
        definition.linear_velocity = {100, 0};
        definition.configs.ccd_override = mode;
        const auto body = mixed.body_pool.add_rigid(definition);
        const auto shape = mixed.shape_pool.add_circle(circle);
        REQUIRE(mixed.body_pool.attach_shape(body, shape, mixed.shape_pool));
        REQUIRE(mixed.island_pool.add_body(island, body, mixed.body_pool));
    }
    REQUIRE(sqp::reschedule_island(mixed, island));
    CHECK(mixed.schedule.entries[island.index].period == short_period);
}

TEST_CASE("CCD islands join existing shorter intervals without setting their deadlines",
          "[world][step][ccd]") {
    sqp::WorldState state;
    state.settings.gravity = {};
    const auto add = [&](sqp::BodyConfigs::CcdOverride mode) {
        sqp::RigidBodyDef definition;
        definition.linear_velocity = {100, 0};
        definition.configs.ccd_override = mode;
        const auto body = state.body_pool.add_rigid(definition);
        sqp::CircleShapeDef circle;
        circle.radius = 1;
        const auto shape = state.shape_pool.add_circle(circle);
        REQUIRE(state.body_pool.attach_shape(body, shape, state.shape_pool));
        const auto island = state.island_pool.create();
        REQUIRE(state.island_pool.add_body(island, body, state.body_pool));
        REQUIRE(sqp::reschedule_island(state, island));
        return island;
    };
    const auto ccd = add(sqp::BodyConfigs::CcdOverride::Auto);
    const auto ordinary = add(sqp::BodyConfigs::CcdOverride::Disabled);
    CHECK(state.schedule.entries[ccd.index].period == sqp::config::TU_PER_TICK);
    CHECK(state.schedule.entries[ordinary.index].period == state.settings.interval.min_period_tu);
    REQUIRE(sqp::run_collect_tasks(state));
    CHECK(state.interval.target_end == state.settings.interval.min_period_tu);
    CHECK(std::find(state.interval.islands->begin(), state.interval.islands->end(), ccd) !=
          state.interval.islands->end());
    CHECK(std::find(state.interval.islands->begin(), state.interval.islands->end(), ordinary) !=
          state.interval.islands->end());
}

TEST_CASE("Collect tasks visits only due period buckets and skips sleep", "[world][step]") {
    sqp::WorldState state;
    state.settings.gravity = {};
    const auto make = [&](sqp::real velocity) {
        sqp::RigidBodyDef definition;
        definition.linear_velocity = {velocity, 0};
        const auto body = state.body_pool.add_rigid(definition);
        const auto id = state.island_pool.create();
        REQUIRE(state.island_pool.add_body(id, body, state.body_pool));
        REQUIRE(sqp::reschedule_island(state, id));
        return id;
    };
    const auto quarter = make(100);
    const auto half = make(10);
    const auto tick = make(1);
    const auto sleeping = make(0);
    REQUIRE(sqp::park_island(state, sleeping));
    const auto bucket_of = [&](sqp::IslandId id) {
        sqp::u8 bucket = 255;
        std::size_t position = 0;
        REQUIRE(state.island_pool.get_schedule_location(id, bucket, position));
        return bucket;
    };
    CHECK(bucket_of(sleeping) == 0);
    CHECK(bucket_of(quarter) < bucket_of(half));
    CHECK(bucket_of(half) < bucket_of(tick));
    const auto due = [&](sqp::tick_t time) {
        state.schedule.time = time - state.settings.interval.min_period_tu;
        REQUIRE(sqp::run_collect_tasks(state));
        CHECK(state.interval.target_end == time);
        CHECK(std::find(state.interval.islands->begin(), state.interval.islands->end(), sleeping) ==
              state.interval.islands->end());
    };
    due(sqp::config::TU_PER_TICK / 4);
    CHECK(state.interval.islands->size() == 1);
    due(sqp::config::TU_PER_TICK / 2);
    CHECK(state.interval.islands->size() == 2);
    due(sqp::config::TU_PER_TICK);
    CHECK(state.interval.islands->size() == 3);
}

TEST_CASE("A changed contact rebuilds only its affected islands", "[world][island]") {
    sqp::WorldState state;
    state.settings.gravity = {};
    std::vector<sqp::IslandId> sleepers;
    for (int i = 0; i < 100; ++i) {
        sqp::RigidBodyDef definition;
        definition.position = {sqp::real(10 + 3 * i), 0};
        const auto body = state.body_pool.add_rigid(definition);
        const auto island = state.island_pool.create();
        REQUIRE(state.island_pool.add_body(island, body, state.body_pool));
        REQUIRE(sqp::park_island(state, island));
        sleepers.push_back(island);
    }
    sqp::CircleShapeDef circle;
    circle.radius = 1;
    const auto make = [&](sqp::real x) {
        sqp::RigidBodyDef definition;
        definition.position = {x, 0};
        const auto body = state.body_pool.add_rigid(definition);
        const auto shape = state.shape_pool.add_circle(circle);
        REQUIRE(state.body_pool.attach_shape(body, shape, state.shape_pool));
        const auto island = state.island_pool.create();
        REQUIRE(state.island_pool.add_body(island, body, state.body_pool));
        REQUIRE(sqp::reschedule_island(state, island));
        return std::pair{body, shape};
    };
    const auto [a, shape_a] = make(0);
    const auto [b, shape_b] = make(1);
    REQUIRE(sqp::run_rebuild_islands(state, false));
    state.contact_graph.set_contact({a, b, shape_a, shape_b, 0, 0}, true);
    REQUIRE(sqp::run_rebuild_islands(state, false));
    CHECK(state.connectivity.last_rebuilt_body_count == 2);
    CHECK(state.body_pool.island_of(a) == state.body_pool.island_of(b));
    for (auto id : sleepers)
        CHECK(state.island_pool.is_sleeping(id));
}

TEST_CASE("Incremental rebuild gives split components distinct island IDs", "[world][island]") {
    sqp::WorldState state;
    const auto a = state.body_pool.add_rigid();
    const auto b = state.body_pool.add_rigid();
    const auto island = state.island_pool.create();
    REQUIRE(state.island_pool.add_body(island, a, state.body_pool));
    REQUIRE(state.island_pool.add_body(island, b, state.body_pool));
    REQUIRE(sqp::reschedule_island(state, island));

    state.connectivity.full_rebuild = false;
    state.connectivity.dirty_islands.push_back(island);
    REQUIRE(sqp::run_rebuild_islands(state, false));

    CHECK(state.body_pool.island_of(a) == island);
    const auto other = state.body_pool.island_of(b);
    REQUIRE(other);
    CHECK(*other != island);
    CHECK(state.island_pool.get(*other) != nullptr);
}

TEST_CASE("Removing one member keeps the remaining island scheduled", "[world][island]") {
    sqp::World world;
    auto& state = world.state();
    const auto a = state.body_pool.add_rigid();
    const auto b = state.body_pool.add_rigid();
    const auto island = state.island_pool.create();
    REQUIRE(state.island_pool.add_body(island, a, state.body_pool));
    REQUIRE(state.island_pool.add_body(island, b, state.body_pool));
    REQUIRE(sqp::reschedule_island(state, island));
    const auto removal = world.destroy_body(0, a);
    REQUIRE(world.step(0));
    CHECK(take(world, removal).status == sqp::CommandStatus::Applied);
    const auto remaining = state.body_pool.island_of(b);
    REQUIRE(remaining);
    sqp::u8 bucket = 255;
    std::size_t position = 0;
    REQUIRE(state.island_pool.get_schedule_location(*remaining, bucket, position));
    CHECK(bucket != 255);
}

TEST_CASE("Persistent contact edges merge and split islands", "[world][island]") {
    sqp::WorldState state;
    sqp::RigidBodyDef definition;
    const auto a = state.body_pool.add_rigid(definition);
    const auto b = state.body_pool.add_rigid(definition);
    sqp::CircleShapeDef circle;
    circle.radius = 1;
    const auto shape_a = state.shape_pool.add_circle(circle);
    const auto shape_b = state.shape_pool.add_circle(circle);
    const auto shape_a2 = state.shape_pool.add_circle(circle);
    const auto shape_b2 = state.shape_pool.add_circle(circle);
    REQUIRE(state.body_pool.attach_shape(a, shape_a, state.shape_pool));
    REQUIRE(state.body_pool.attach_shape(b, shape_b, state.shape_pool));
    REQUIRE(state.body_pool.attach_shape(a, shape_a2, state.shape_pool));
    REQUIRE(state.body_pool.attach_shape(b, shape_b2, state.shape_pool));
    const auto island_a = state.island_pool.create();
    const auto island_b = state.island_pool.create();
    REQUIRE(state.island_pool.add_body(island_a, a, state.body_pool));
    REQUIRE(state.island_pool.add_body(island_b, b, state.body_pool));
    REQUIRE(state.island_pool.set_toi_truncations(island_a, 2));
    REQUIRE(state.island_pool.set_toi_truncations(island_b, 3));
    state.contact_graph.set_contact({a, b, shape_a, shape_b, 0, 0}, true);
    state.contact_graph.set_contact({a, b, shape_a2, shape_b2, 0, 0}, true);
    REQUIRE(sqp::run_rebuild_islands(state));
    CHECK(state.body_pool.island_of(a) == state.body_pool.island_of(b));
    CHECK(state.body_pool.island_of(a) == island_a);
    CHECK(state.island_pool.toi_truncations(island_a) == 3);

    REQUIRE(state.body_pool.detach_shape(a, shape_a, state.shape_pool));
    REQUIRE(state.shape_pool.remove(shape_a));
    REQUIRE(sqp::run_rebuild_islands(state));
    CHECK(state.contact_graph.contacts().size() == 1);
    CHECK(state.body_pool.island_of(a) == state.body_pool.island_of(b));

    REQUIRE(state.body_pool.set_pose(b, {10, 0}, 0));
    REQUIRE(state.island_pool.set_toi_truncations(island_a, 4));
    REQUIRE(sqp::run_rebuild_islands(state));
    CHECK(state.contact_graph.contacts().empty());
    CHECK(state.body_pool.island_of(a) != state.body_pool.island_of(b));
    CHECK(state.body_pool.island_of(a) == island_a);
    CHECK(state.island_pool.toi_truncations(*state.body_pool.island_of(a)) == 4);
    CHECK(state.island_pool.toi_truncations(*state.body_pool.island_of(b)) == 4);
}

TEST_CASE("Contact adjacency survives edge swaps and body slot reuse", "[world][island]") {
    sqp::ConstraintGraph graph;
    using Kind = sqp::BodyConfigs::Kind;
    const sqp::BodyHandle a{0, 1, Kind::Rigid};
    const sqp::BodyHandle b{1, 1, Kind::Rigid};
    const sqp::BodyHandle c{2, 1, Kind::Rigid};
    const sqp::BodyHandle replacement{1, 2, Kind::Rigid};
    const sqp::ShapeHandle sa{0, 1, sqp::ShapeKind::Circle};
    const sqp::ShapeHandle sb{1, 1, sqp::ShapeKind::Circle};
    const sqp::ShapeHandle sc{2, 1, sqp::ShapeKind::Circle};
    const sqp::ShapeHandle sr{1, 2, sqp::ShapeKind::Circle};
    const sqp::Pair ab{a, b, sa, sb, 0, 0};
    const sqp::Pair ac{a, c, sa, sc, 0, 0};
    const sqp::Pair ar{a, replacement, sa, sr, 0, 0};

    graph.set_contact(ab, true);
    graph.set_contact(ac, true);
    graph.set_contact(ab, false); // Moves ac into ab's former contact index.
    REQUIRE(graph.contacts().size() == 1);
    CHECK(graph.contacts_for(c).size() == 1);
    graph.set_contact(ac, false);
    CHECK(graph.contacts().empty());

    graph.set_contact(ab, true);
    graph.set_contact(ar, true); // Reuses b's slot with a new generation.
    CHECK(graph.contacts_for(b).empty());
    CHECK(graph.contacts_for(replacement).size() == 1);
    graph.erase_if([&](const auto& edge) { return edge.body_b == b; });
    REQUIRE(graph.contacts().size() == 1);
    CHECK(graph.contacts_for(replacement).size() == 1);
    graph.set_contact(ar, false);
    CHECK(graph.contacts().empty());
}

TEST_CASE("Contact connectivity joins rigid, lite circle and AABB bodies", "[world][island]") {
    sqp::WorldState state;
    const auto rigid = state.body_pool.add_rigid();
    const auto rigid_shape = state.shape_pool.add_circle({1});
    const auto circle_shape = state.shape_pool.add_circle({1});
    sqp::PolygonShapeDef box_def;
    box_def.vertices = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
    const auto box_shape = state.shape_pool.add_polygon(box_def);
    REQUIRE(state.body_pool.attach_shape(rigid, rigid_shape, state.shape_pool));
    const auto circle = state.body_pool.add_lite_circle({}, circle_shape, state.shape_pool);
    const auto box = state.body_pool.add_aabb_box({}, box_shape, state.shape_pool);
    REQUIRE(state.body_pool.is_valid(rigid));
    REQUIRE(state.body_pool.is_valid(circle));
    REQUIRE(state.body_pool.is_valid(box));
    const auto rigid_island = state.island_pool.create();
    const auto circle_island = state.island_pool.create();
    const auto box_island = state.island_pool.create();
    REQUIRE(state.island_pool.add_body(rigid_island, rigid, state.body_pool));
    REQUIRE(state.island_pool.add_body(circle_island, circle, state.body_pool));
    REQUIRE(state.island_pool.add_body(box_island, box, state.body_pool));
    state.contact_graph.set_contact({rigid, circle, rigid_shape, circle_shape, 0, 0}, true);
    state.contact_graph.set_contact({circle, box, circle_shape, box_shape, 0, 0}, true);

    REQUIRE(sqp::run_rebuild_islands(state));
    CHECK(state.contact_graph.contacts().size() == 2);
    CHECK(state.body_pool.island_of(rigid) == rigid_island);
    CHECK(state.body_pool.island_of(circle) == rigid_island);
    CHECK(state.body_pool.island_of(box) == rigid_island);
    REQUIRE(sqp::run_rebuild_islands(state));
    CHECK(state.body_pool.island_of(box) == rigid_island);
}

TEST_CASE("Connectivity lookup follows dense moves and rejects reused body slots",
          "[world][island]") {
    sqp::WorldState state;
    const auto removed = state.body_pool.add_rigid();
    const auto a = state.body_pool.add_rigid();
    const auto b = state.body_pool.add_rigid();
    const auto shape_a = state.shape_pool.add_circle({1});
    const auto shape_b = state.shape_pool.add_circle({1});
    REQUIRE(state.body_pool.attach_shape(a, shape_a, state.shape_pool));
    REQUIRE(state.body_pool.attach_shape(b, shape_b, state.shape_pool));
    const auto removed_island = state.island_pool.create();
    const auto island_a = state.island_pool.create();
    const auto island_b = state.island_pool.create();
    REQUIRE(state.island_pool.add_body(removed_island, removed, state.body_pool));
    REQUIRE(state.island_pool.add_body(island_a, a, state.body_pool));
    REQUIRE(state.island_pool.add_body(island_b, b, state.body_pool));
    state.contact_graph.set_contact({a, b, shape_a, shape_b, 0, 0}, true);

    REQUIRE(state.body_pool.remove(removed, state.island_pool));
    REQUIRE(state.body_pool.dense_index(b));
    CHECK(*state.body_pool.dense_index(b) == 0);
    REQUIRE(sqp::run_rebuild_islands(state));
    CHECK(state.body_pool.island_of(a) == state.body_pool.island_of(b));

    REQUIRE(state.body_pool.remove(a, state.island_pool));
    const auto replacement = state.body_pool.add_rigid();
    CHECK(replacement.index() == a.index());
    CHECK(replacement.generation != a.generation);
    const auto replacement_shape = state.shape_pool.add_circle({1});
    REQUIRE(state.body_pool.attach_shape(replacement, replacement_shape, state.shape_pool));
    const auto replacement_island = state.island_pool.create();
    REQUIRE(state.island_pool.add_body(replacement_island, replacement, state.body_pool));
    REQUIRE(sqp::run_rebuild_islands(state));
    CHECK(state.contact_graph.contacts().empty());
    CHECK(state.body_pool.island_of(replacement) != state.body_pool.island_of(b));

    state.contact_graph.set_contact({replacement, b, replacement_shape, shape_b, 0, 0}, true);
    REQUIRE(sqp::run_rebuild_islands(state));
    CHECK(state.body_pool.island_of(replacement) == state.body_pool.island_of(b));
    REQUIRE(sqp::run_rebuild_islands(state));
    CHECK(state.body_pool.island_of(replacement) == state.body_pool.island_of(b));
}

TEST_CASE("Tick interval resets TOI budgets after shorter ordinary intervals", "[world][step]") {
    sqp::WorldState state;
    CHECK(state.schedule.period_islands.size() == sqp::config::PERIOD_BUCKET_COUNT + 1);
    state.settings.gravity = {};
    sqp::RigidBodyDef definition;
    definition.linear_velocity = {100, 0};
    const auto body = state.body_pool.add_rigid(definition);
    const auto island = state.island_pool.create();
    REQUIRE(state.island_pool.add_body(island, body, state.body_pool));
    REQUIRE(state.island_pool.set_toi_truncations(island, 2));
    const sqp::tick_t period = state.settings.interval.min_period_tu;
    for (sqp::tick_t time = period; time <= sqp::config::TU_PER_TICK; time += period) {
        REQUIRE(sqp::run_collect_tasks(state));
        CHECK(state.interval.target_end == time);
        CHECK(state.interval.tick == (time == sqp::config::TU_PER_TICK));
        REQUIRE(sqp::run_advance_interval(state));
        CHECK(state.island_pool.toi_truncations(island) ==
              (time == sqp::config::TU_PER_TICK ? 0 : 2));
    }
}

TEST_CASE("World stage statistics count dispatched calls and completed intervals",
          "[world][profile]") {
    using Stage = sqp::WorldState::StageStats::Stage;
    const auto index = [](Stage stage) { return static_cast<std::size_t>(stage); };
    sqp::World world;
    const auto creation = world.create_rigid(0, {});
    REQUIRE(world.step(0));
    REQUIRE(take(world, creation).status == sqp::CommandStatus::Applied);
    world.state().stats = {};
    REQUIRE(world.step(sqp::config::TU_PER_TICK));
    const auto& stats = world.state().stats;
    CHECK(stats.intervals > 0);
    CHECK(stats.stages[index(Stage::CollectTasks)].calls > 0);
    CHECK(stats.stages[index(Stage::AdvanceInterval)].calls == stats.intervals);
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
    CHECK(stats.stages[index(Stage::CollectTasks)].nanoseconds > 0);
    REQUIRE(stats.interval_entries.size() == stats.intervals);
    CHECK(stats.interval_entries.front().active_islands == 1);
    CHECK(stats.interval_entries.front().active_bodies == 1);
#else
    CHECK(stats.stages[index(Stage::CollectTasks)].nanoseconds == 0);
    CHECK(stats.stages[index(Stage::AdvanceInterval)].nanoseconds == 0);
#endif

    world.state().stats = {};
    CHECK(world.state().stats.intervals == 0);
    CHECK(world.state().stats.stages[index(Stage::AdvanceInterval)].calls == 0);
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
    CHECK(world.state().stats.interval_entries.empty());
#endif
}

TEST_CASE("A static contact endpoint never joins two dynamic islands", "[world][island]") {
    sqp::WorldState state;
    sqp::RigidBodyDef dynamic, fixed;
    fixed.configs.register_type = sqp::BodyConfigs::RegisterType::ImmutableStatic;
    const auto a = state.body_pool.add_rigid(dynamic);
    const auto b = state.body_pool.add_rigid(dynamic);
    const auto ground = state.body_pool.add_rigid(fixed);
    sqp::CircleShapeDef circle;
    circle.radius = 1;
    const auto shape_a = state.shape_pool.add_circle(circle);
    const auto shape_b = state.shape_pool.add_circle(circle);
    const auto shape_ground = state.shape_pool.add_circle(circle);
    REQUIRE(state.body_pool.attach_shape(a, shape_a, state.shape_pool));
    REQUIRE(state.body_pool.attach_shape(b, shape_b, state.shape_pool));
    REQUIRE(state.body_pool.attach_shape(ground, shape_ground, state.shape_pool));
    const auto island_a = state.island_pool.create();
    const auto island_b = state.island_pool.create();
    REQUIRE(state.island_pool.add_body(island_a, a, state.body_pool));
    REQUIRE(state.island_pool.add_body(island_b, b, state.body_pool));
    state.contact_graph.set_contact({a, ground, shape_a, shape_ground, 0, 0}, true);
    state.contact_graph.set_contact({b, ground, shape_b, shape_ground, 0, 0}, true);
    REQUIRE(sqp::run_rebuild_islands(state));
    CHECK(state.body_pool.island_of(a) != state.body_pool.island_of(b));
    CHECK_FALSE(state.body_pool.island_of(ground));
}

TEST_CASE("World consumes a TOI microstep and preserves later ordinary work",
          "[world][step][toi]") {
    sqp::World world;
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    world.set_settings(0, settings);
    sqp::RigidBodyDef moving;
    moving.position = {0, 0};
    moving.linear_velocity = {1000, 0};
    moving.configs.event_flags = sqp::BodyConfigs::EventFlags::AfterSolve;
    sqp::RigidBodyDef stationary;
    stationary.position = {4, 0};
    const auto moving_id = world.create_rigid(0, moving);
    const auto stationary_id = world.create_rigid(0, stationary);
    int events = 0;
    REQUIRE(world.step(0));
    const auto a = take(world, moving_id).body;
    const auto b = take(world, stationary_id).body;
    sqp::CircleShapeDef circle;
    circle.radius = sqp::real(0.5);
    const auto shape_a = world.create_circle_shape(0, circle, a);
    const auto shape_b = world.create_circle_shape(0, circle, b);
    REQUIRE(world.step(0));
    REQUIRE(take(world, shape_a).status == sqp::CommandStatus::Applied);
    REQUIRE(take(world, shape_b).status == sqp::CommandStatus::Applied);
    REQUIRE_FALSE(world.state().process_toi.toi_heap.empty());
    const auto toi_time = world.state().process_toi.toi_heap.front().time;
    REQUIRE(toi_time > 20000);
    REQUIRE(toi_time < sqp::config::TU_PER_TICK);
    const auto registration = world.set_after_solve_callback(toi_time, on_after_solve, &events);
    REQUIRE(world.step(20000));
    CHECK(body_state(world, a).position.x == 0);
    CHECK(world.state().pairs.delayed_pairs->empty());
    const auto before_microstep = world.state().stats.intervals;
    REQUIRE(world.step(toi_time - 20000));
    CHECK(world.state().stats.intervals > before_microstep);
    CHECK(world.time_tu() == toi_time);
    const auto before_tick_island = world.query().bodies().island_of(a);
    REQUIRE(before_tick_island);
    CHECK(world.state().island_pool.toi_truncations(*before_tick_island) >= 1);
    const auto after_microstep = world.state().stats.intervals;
    REQUIRE(world.step(sqp::config::TU_PER_TICK - toi_time));
    CHECK(world.state().stats.intervals > after_microstep);
    CHECK(take(world, registration).status == sqp::CommandStatus::Applied);
    CHECK(world.time_tu() == sqp::config::TU_PER_TICK);
    CHECK_FALSE(world.state().interval.microstep);
    CHECK(world.state().pairs.delayed_pairs->empty());
    CHECK(world.state().interval.delayed_islands->empty());
    CHECK(events >= 1);
    const auto impact_island = world.query().bodies().island_of(a);
    REQUIRE(impact_island);
    CHECK(world.state().island_pool.toi_truncations(*impact_island) == 0);
    REQUIRE(world.step(sqp::config::TU_PER_TICK));
    CHECK(world.time_tu() == sqp::config::TU_PER_TICK * 2);
}

TEST_CASE("A cancelled TOI leaves the body its full tick of motion", "[world][step][toi]") {
    using namespace sqp;
    World world;
    PhysicsSettings settings;
    settings.gravity = {};
    settings.interval.max_truncations_per_tick = 1;
    const auto settings_id = world.set_settings(0, settings);
    RigidBodyDef moving({0, 0}, 0, {1024, 0}, 0);
    moving.configs.event_flags = BodyConfigs::EventFlags::PreSolve;
    RigidBodyDef fixed;
    fixed.position = {4, 0};
    fixed.configs.register_type = BodyConfigs::RegisterType::ImmutableStatic;
    const auto moving_id = world.create_rigid(0, moving);
    const auto fixed_id = world.create_rigid(0, fixed);

    int cancelled = 0;
    const auto cancel = [](PhysicsEvent& event, void* data) {
        auto& contact = std::get<PreSolveEvent>(event.event);
        contact.manifold.count = 0;
        ++*static_cast<int*>(data);
    };
    const auto callback_id = world.set_pre_solve_callback(0, cancel, &cancelled);
    REQUIRE(world.step(0));
    REQUIRE(take(world, settings_id).status == CommandStatus::Applied);
    REQUIRE(take(world, callback_id).status == CommandStatus::Applied);
    const auto body = take(world, moving_id).body;
    const auto obstacle = take(world, fixed_id).body;
    CircleShapeDef circle;
    circle.radius = real(0.5);
    const auto shape_a = world.create_circle_shape(0, circle, body);
    const auto shape_b = world.create_circle_shape(0, circle, obstacle);
    REQUIRE(world.step(config::TU_PER_TICK / 2));
    REQUIRE(take(world, shape_a).status == CommandStatus::Applied);
    REQUIRE(take(world, shape_b).status == CommandStatus::Applied);
    REQUIRE(cancelled > 0);
    const auto island = world.query().bodies().island_of(body);
    REQUIRE(island);
    CHECK(world.state().island_pool.toi_truncations(*island) == 0);

    REQUIRE(world.step(config::TU_PER_TICK / 2));
    CHECK(body_state(world, body).position.x == Catch::Approx(16));
    CHECK(world.state().island_pool.get(*island)->time == config::TU_PER_TICK);
    CHECK(world.state().island_pool.last_frozen_until(*island) == 0);
}

TEST_CASE("Static bodies keep proxies without owning an island", "[world][step]") {
    sqp::World world;
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    world.set_settings(0, settings);
    sqp::RigidBodyDef fixed, moving;
    fixed.configs.register_type = sqp::BodyConfigs::RegisterType::ImmutableStatic;
    fixed.position = {0, 0};
    moving.position = {0, sqp::real(1.5)};
    const auto fixed_id = world.create_rigid(0, fixed);
    const auto moving_id = world.create_rigid(0, moving);
    REQUIRE(world.step(0));
    const auto ground = take(world, fixed_id).body;
    const auto body = take(world, moving_id).body;
    CHECK_FALSE(world.query().bodies().island_of(ground));
    REQUIRE(world.query().bodies().island_of(body));
    sqp::CircleShapeDef circle;
    circle.radius = 1;
    world.create_circle_shape(0, circle, ground);
    world.create_circle_shape(0, circle, body);
    REQUIRE(world.step(0));
    REQUIRE(world.step(sqp::config::TU_PER_TICK));
    CHECK_FALSE(world.query().bodies().island_of(ground));
}

TEST_CASE("A quiet island sleeps and a timed velocity command wakes it", "[world][step]") {
    sqp::World world;
    sqp::PhysicsSettings settings;
    settings.gravity = {};
    settings.sleep_time_threshold_tu = sqp::config::TU_PER_TICK;
    world.set_settings(0, settings);
    const auto creation = world.create_rigid(0);
    REQUIRE(world.step(0));
    const auto body = take(world, creation).body;
    REQUIRE(world.step(sqp::config::TU_PER_TICK));
    CHECK(world.state().island_pool.is_sleeping(*world.query().bodies().island_of(body)));
    const auto island = world.query().bodies().island_of(body);
    REQUIRE(island);
    CHECK(world.state().island_pool.sleep_time(*island) >= sqp::config::TU_PER_TICK);

    const auto change = world.set_body_velocity(world.time_tu(), body, {1, 0});
    REQUIRE(world.step(0));
    CHECK(take(world, change).status == sqp::CommandStatus::Applied);
    CHECK_FALSE(world.state().island_pool.is_sleeping(*island));
    CHECK(world.state().island_pool.sleep_time(*island) == 0);
}
