#include "world/world_stages.h"

#include <catch_amalgamated.hpp>

#include <limits>
#include <vector>

namespace {

struct Fixture {
    sqp::WorldState state;
    sqp::BodyHandle ground, moving;
    sqp::ShapeHandle ground_shape, moving_shape;
    sqp::IslandId island;

    Fixture(sqp::Vec2 velocity = {0, -2}) {
        state.settings.gravity = {};
        sqp::RigidBodyDef ground_def;
        ground_def.configs.register_type = sqp::BodyConfigs::RegisterType::ImmutableStatic;
        ground = state.body_pool.add_rigid(ground_def);
        moving = state.body_pool.add_rigid({{0, 1}, 0, velocity, 0});
        ground_shape = state.shape_pool.add_circle({sqp::real(0.5)});
        moving_shape = state.shape_pool.add_circle({sqp::real(0.5)});
        REQUIRE(state.body_pool.attach_shape(ground, ground_shape, state.shape_pool));
        REQUIRE(state.body_pool.attach_shape(moving, moving_shape, state.shape_pool));
        island = state.island_pool.create();
        REQUIRE(state.island_pool.add_body(island, ground, state.body_pool));
        REQUIRE(state.island_pool.add_body(island, moving, state.body_pool));
        state.interval.islands->push_back(island);
        state.interval.target_end = sqp::config::TU_PER_SECOND / 4;
    }

    void append(sqp::real depth = sqp::real(0.01)) {
        const auto ia = *state.body_pool.dense_index(ground);
        const auto ib = *state.body_pool.dense_index(moving);
        state.pairs.pairs->push_back({ground, moving, ground_shape, moving_shape,
                                      state.body_pool.motion_version_column(ground.kind())[ia],
                                      state.body_pool.motion_version_column(moving.kind())[ib]});
        sqp::ContactManifold manifold;
        manifold.count = 1;
        manifold.normal = {0, 1};
        manifold.points[0].point = {0, sqp::real(0.5)};
        manifold.points[0].depth = depth;
        REQUIRE(state.pairs.manifolds->append(manifold));
    }

    void solve() {
        REQUIRE(sqp::run_pre_solve(state));
        REQUIRE(sqp::run_solve(state));
    }
};

struct CallbackData {
    std::vector<sqp::PostSolveEvent> received;
    int pre_calls = 0;
    bool disable_first = false;
};

void pre_callback(sqp::PhysicsEvent& event, void* context) {
    auto& data = *static_cast<CallbackData*>(context);
    auto& contact = std::get<sqp::PreSolveEvent>(event.event);
    ++data.pre_calls;
    contact.manifold.points[0].point.x = sqp::real(0.25);
    if (data.disable_first && data.pre_calls == 1)
        contact.manifold.count = 0;
}

void post_callback(sqp::PhysicsEvent& event, void* context) {
    auto& data = *static_cast<CallbackData*>(context);
    REQUIRE(event.type == sqp::PhysicsEventTypes::PostSolve);
    auto& contact = std::get<sqp::PostSolveEvent>(event.event);
    data.received.push_back(contact);
    contact.manifold.points[0].point.x = 99;
    contact.report.total_impulse = 99;
}

} // namespace

TEST_CASE("PostSolve reports original manifolds in solved pair order and clears consumed data",
          "[world][post_solve]") {
    Fixture f;
    f.append();
    f.append();
    CallbackData data;
    f.state.callbacks = {pre_callback, post_callback, &data, &data};
    REQUIRE(
        f.state.body_pool.set_event_flags(f.moving, sqp::BodyConfigs::EventFlags::PreSolve |
                                                        sqp::BodyConfigs::EventFlags::AfterSolve));
    f.solve();
    REQUIRE(f.state.post_solve.reports.size() == 2);
    CHECK(f.state.post_solve.reports[0].pair_index == 0);
    CHECK(f.state.post_solve.reports[1].pair_index == 1);
    CHECK(f.state.solve.constraint_pair_indices.empty());
    CHECK(f.state.pairs.combined_restitutions.empty());
    const auto active = f.state.pairs.manifolds->at(0);
    REQUIRE(active.points[0].point.x == sqp::real(0.25));

    REQUIRE(sqp::run_post_solve(f.state));
    REQUIRE(data.received.size() == 2);
    for (std::size_t i = 0; i < data.received.size(); ++i) {
        CHECK(data.received[i].contact_pair == (*f.state.pairs.pairs)[i]);
        CHECK(data.received[i].manifold.points[0].point.x == 0);
        CHECK(data.received[i].report.normal_impulse >= 0);
    }
    CHECK(f.state.pairs.manifolds->at(0).points[0].point.x == sqp::real(0.25));
    CHECK(f.state.pairs.pairs->size() == 2);
    CHECK_FALSE(f.state.pairs.original_manifolds_valid);
    CHECK(f.state.pairs.combined_restitutions.empty());
    CHECK(f.state.pairs.combined_contact_frequencies.empty());
    CHECK(f.state.solve.constraints.empty());
    CHECK(f.state.solve.prepared.empty());
    CHECK(f.state.solve.constraint_pair_indices.empty());
    CHECK(f.state.post_solve.reports.empty());
    CHECK(f.state.post_solve.events.empty());
}

TEST_CASE("PostSolve skips disabled contacts but reports zero impulse constraints",
          "[world][post_solve]") {
    SECTION("disabled first row") {
        Fixture f;
        f.append();
        f.append();
        CallbackData data;
        data.disable_first = true;
        f.state.callbacks = {pre_callback, post_callback, &data, &data};
        REQUIRE(f.state.body_pool.set_event_flags(f.ground,
                                                  sqp::BodyConfigs::EventFlags::PreSolve |
                                                      sqp::BodyConfigs::EventFlags::AfterSolve));
        f.solve();
        REQUIRE(f.state.post_solve.reports.size() == 1);
        CHECK(f.state.post_solve.reports[0].pair_index == 1);
        REQUIRE(sqp::run_post_solve(f.state));
        REQUIRE(data.received.size() == 1);
        CHECK(data.received[0].contact_pair == (*f.state.pairs.pairs)[1]);
    }
    SECTION("zero impulse still reports") {
        Fixture f({0, 0});
        f.append(0);
        CallbackData data;
        f.state.callbacks.post_solve = post_callback;
        f.state.callbacks.post_data = &data;
        REQUIRE(
            f.state.body_pool.set_event_flags(f.moving, sqp::BodyConfigs::EventFlags::AfterSolve));
        f.solve();
        REQUIRE(f.state.post_solve.reports.size() == 1);
        CHECK(f.state.post_solve.reports[0].pair_index == 0);
        CHECK(f.state.post_solve.reports[0].contact.normal_impulse == 0);
        REQUIRE(sqp::run_post_solve(f.state));
        REQUIRE(data.received.size() == 1);
        CHECK(data.received[0].report.normal_impulse == 0);
    }
}

TEST_CASE("PostSolve respects callback registration and body flags", "[world][post_solve]") {
    for (int mode = 0; mode < 5; ++mode) {
        Fixture f;
        f.append();
        CallbackData data;
        if (mode != 0 && mode != 4) {
            f.state.callbacks.post_solve = post_callback;
            f.state.callbacks.post_data = &data;
        }
        if (mode >= 2)
            REQUIRE(f.state.body_pool.set_event_flags(f.ground,
                                                      sqp::BodyConfigs::EventFlags::AfterSolve));
        if (mode == 3)
            REQUIRE(f.state.body_pool.set_event_flags(f.moving,
                                                      sqp::BodyConfigs::EventFlags::AfterSolve));
        f.solve();
        CHECK(f.state.post_solve.reports.size() == (mode >= 2 ? 1 : 0));
        REQUIRE(sqp::run_post_solve(f.state));
        CHECK(data.received.size() == (mode == 2 || mode == 3 ? 1 : 0));
    }
}

TEST_CASE("PostSolve reports and cleans ordinary intervals and microsteps", "[world][post_solve]") {
    for (bool microstep : {false, true}) {
        Fixture f;
        f.append();
        CallbackData data;
        f.state.callbacks = {pre_callback, post_callback, &data, &data};
        REQUIRE(f.state.body_pool.set_event_flags(f.moving,
                                                  sqp::BodyConfigs::EventFlags::PreSolve |
                                                      sqp::BodyConfigs::EventFlags::AfterSolve));
        f.state.interval.microstep = microstep;
        f.state.interval.microstep_time = sqp::config::TU_PER_SECOND / 2;
        f.solve();
        REQUIRE(sqp::run_post_solve(f.state));
        REQUIRE(data.received.size() == 1);
        CHECK(data.received[0].manifold.points[0].point.x == 0);
        CHECK_FALSE(f.state.pairs.original_manifolds_valid);
        CHECK(f.state.post_solve.reports.empty());
        CHECK(f.state.pairs.combined_restitutions.empty());
    }
}

#ifdef SQP_DEBUG_VALIDATION
TEST_CASE("PostSolve validates all rows before sending any callback", "[world][post_solve]") {
    Fixture f;
    f.append();
    f.append();
    CallbackData data;
    f.state.callbacks.post_solve = post_callback;
    f.state.callbacks.post_data = &data;
    REQUIRE(f.state.body_pool.set_event_flags(f.moving, sqp::BodyConfigs::EventFlags::AfterSolve));
    f.solve();
    REQUIRE(f.state.post_solve.reports.size() == 2);
    const auto original = f.state.post_solve.reports[1].contact.normal_impulse;
    f.state.post_solve.reports[1].contact.normal_impulse =
        std::numeric_limits<sqp::real>::quiet_NaN();
    CHECK_FALSE(sqp::run_post_solve(f.state));
    CHECK(data.received.empty());
    CHECK(f.state.pairs.original_manifolds_valid);
    CHECK(f.state.post_solve.reports.size() == 2);
    CHECK(f.state.pairs.combined_restitutions.empty());
    f.state.post_solve.reports[1].contact.normal_impulse = original;
    f.state.post_solve.reports[1].pair_index = 0;
    CHECK_FALSE(sqp::run_post_solve(f.state));
    CHECK(data.received.empty());
    f.state.post_solve.reports[1].pair_index = 1;
    const auto shape = f.state.pairs.pairs->at(1).shape_b;
    f.state.pairs.pairs->at(1).shape_b = f.ground_shape;
    CHECK_FALSE(sqp::run_post_solve(f.state));
    CHECK(data.received.empty());
    f.state.pairs.pairs->at(1).shape_b = shape;
    REQUIRE(sqp::run_post_solve(f.state));
    CHECK(data.received.size() == 2);
}
#endif

TEST_CASE("PostSolve rejects an out-of-range report before callbacks",
          "[world][post_solve]") {
    Fixture f;
    f.append();
    CallbackData data;
    f.state.callbacks.post_solve = post_callback;
    f.state.callbacks.post_data = &data;
    REQUIRE(f.state.body_pool.set_event_flags(f.moving,
                                               sqp::BodyConfigs::EventFlags::AfterSolve));
    f.solve();
    REQUIRE(f.state.post_solve.reports.size() == 1);
    // Report indices remain checked even without debug validation.
    f.state.post_solve.reports[0].pair_index = f.state.pairs.pairs->size();
    CHECK_FALSE(sqp::run_post_solve(f.state));
    CHECK(data.received.empty());
    CHECK(f.state.pairs.original_manifolds_valid);
}

TEST_CASE("TOI switching parks the original manifold with ordinary pairs",
          "[world][post_solve][toi]") {
    Fixture f;
    f.append();
    REQUIRE(sqp::run_pre_solve(f.state));
    REQUIRE(f.state.pairs.original_manifolds_valid);
    f.state.narrow_phase.toi_pairs.push_back(f.state.pairs.pairs->at(0));
    f.state.narrow_phase.toi_times.push_back(f.state.interval.target_end);
    f.state.narrow_phase.toi_normals.push_back({0, -1});
    bool selected = false;
    REQUIRE(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
    REQUIRE(selected);
    CHECK(f.state.interval.microstep);
    REQUIRE(f.state.pairs.delayed_original_manifolds_valid);
    CHECK(f.state.pairs.delayed_original_manifolds.size() == 1);
    CHECK(f.state.pairs.delayed_original_manifolds.at(0).points[0].point.x == 0);
    CHECK_FALSE(f.state.pairs.original_manifolds_valid);
}

TEST_CASE("TOI selection parks the original manifold when another event is stale",
          "[world][post_solve][toi]") {
    Fixture f;
    f.append();
    REQUIRE(sqp::run_pre_solve(f.state));
    const auto ordinary_pair = f.state.pairs.pairs->at(0);
    const auto late = f.state.body_pool.add_rigid();
    const auto late_shape = f.state.shape_pool.add_circle({sqp::real(0.5)});
    REQUIRE(f.state.body_pool.attach_shape(late, late_shape, f.state.shape_pool));
    const auto late_island = f.state.island_pool.create();
    REQUIRE(f.state.island_pool.add_body(late_island, late, f.state.body_pool));
    REQUIRE(f.state.island_pool.set_time(late_island, f.state.interval.target_end));
    const auto ground_index = *f.state.body_pool.dense_index(f.ground);
    const auto late_index = *f.state.body_pool.dense_index(late);
    const sqp::Pair late_pair{
        f.ground,
        late,
        f.ground_shape,
        late_shape,
        f.state.body_pool.motion_version_column(f.ground.kind())[ground_index],
        f.state.body_pool.motion_version_column(late.kind())[late_index]};
    f.state.narrow_phase.toi_pairs = {ordinary_pair, late_pair};
    f.state.narrow_phase.toi_times = {f.state.interval.target_end, f.state.interval.target_end};
    f.state.narrow_phase.toi_normals = {{0, -1}, {0, -1}};
    bool selected = false;
    CHECK(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
    CHECK(selected);
    CHECK(f.state.pairs.pairs->size() == 1);
    CHECK(f.state.pairs.pairs->at(0) == ordinary_pair);
    CHECK_FALSE(f.state.pairs.original_manifolds_valid);
    CHECK(f.state.pairs.delayed_original_manifolds_valid);
    CHECK(f.state.pairs.delayed_original_manifolds.at(0).points[0].point.x == 0);
}
