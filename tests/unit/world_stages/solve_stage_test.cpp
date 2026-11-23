#include "world/world_stages.h"

#include <catch_amalgamated.hpp>

#include <initializer_list>
#include <limits>
#include <vector>

namespace {

struct Fixture {
    sqp::WorldState state;

    sqp::BodyHandle
    rigid(sqp::Vec2 position = {}, sqp::Vec2 velocity = {},
          sqp::BodyConfigs::RegisterType type = sqp::BodyConfigs::RegisterType::Dynamic,
          sqp::BodyConfigs::OverrideFlags flags = {}) {
        sqp::RigidBodyDef def(position, 0, velocity, 0);
        def.configs.register_type = type;
        def.configs.override_flags = flags;
        const auto body = state.body_pool.add_rigid(def);
        const auto shape = state.shape_pool.add_circle({sqp::real(1)});
        REQUIRE(state.body_pool.attach_shape(body, shape, state.shape_pool));
        return body;
    }

    sqp::IslandId island(std::initializer_list<sqp::BodyHandle> bodies, sqp::tick_t time = 0,
                         bool active = true) {
        const auto id = state.island_pool.create();
        REQUIRE(state.island_pool.set_time(id, time));
        for (auto body : bodies)
            REQUIRE(state.island_pool.add_body(id, body, state.body_pool));
        if (active)
            state.interval.islands->push_back(id);
        return id;
    }

    void contact(sqp::BodyHandle a, sqp::BodyHandle b, sqp::real depth = sqp::real(0.01)) {
        sqp::ShapeHandle shape_a, shape_b;
        std::vector<sqp::ShapeHandle> shapes;
        REQUIRE(state.body_pool.get_attached_shapes(a, shapes));
        REQUIRE(shapes.size() == 1);
        shape_a = shapes[0];
        REQUIRE(state.body_pool.get_attached_shapes(b, shapes));
        REQUIRE(shapes.size() == 1);
        shape_b = shapes[0];
        const auto index_a = *state.body_pool.dense_index(a);
        const auto index_b = *state.body_pool.dense_index(b);
        state.pairs.pairs->push_back({a, b, shape_a, shape_b,
                                      state.body_pool.motion_version_column(a.kind())[index_a],
                                      state.body_pool.motion_version_column(b.kind())[index_b]});
        sqp::ContactManifold manifold;
        manifold.count = 1;
        manifold.normal = {0, 1};
        manifold.points[0].point = {0, sqp::real(0.5)};
        manifold.points[0].depth = depth;
        REQUIRE(state.pairs.manifolds->append(manifold));
    }

    sqp::BodyState body(sqp::BodyHandle handle) const {
        sqp::BodyState result;
        REQUIRE(state.body_pool.get_state(handle, result));
        return result;
    }
};

} // namespace

TEST_CASE("Solve integrates force, torque and gravity without contacts", "[world][solve]") {
    Fixture f;
    f.state.settings.gravity = {0, -10};
    const auto body = f.rigid({}, {1, 1});
    const auto id = f.island({body});
    sqp::BodyProperties properties;
    REQUIRE(f.state.body_pool.get_properties(body, properties));
    REQUIRE(f.state.body_pool.set_constant_force(body, {properties.mass * 2, 0},
                                                 properties.inertia * 4));
    f.state.interval.target_end = sqp::config::TU_PER_SECOND / 2;

    REQUIRE(sqp::run_pre_solve(f.state));
    REQUIRE(sqp::run_solve(f.state));
    const auto after = f.body(body);
    CHECK(after.linear_velocity.x == Catch::Approx(2));
    CHECK(after.linear_velocity.y == Catch::Approx(-4));
    CHECK(after.angular_velocity == Catch::Approx(2));
    CHECK(after.position.x == Catch::Approx(0.8125));
    CHECK(after.position.y == Catch::Approx(-1.0625));
    CHECK(after.rotation == Catch::Approx(0.625));
    CHECK(f.state.island_pool.get(id)->time == f.state.interval.target_end);
    CHECK(f.state.post_solve.reports.empty());
}

TEST_CASE("Ordinary Solve uses floating substep durations", "[world][solve]") {
    for (const sqp::u32 steps : {1u, 4u}) {
        Fixture f;
        f.state.settings.gravity = {};
        f.state.settings.contact.substeps = steps;
        const auto body = f.rigid();
        const auto id = f.island({body});
        sqp::BodyProperties properties;
        REQUIRE(f.state.body_pool.get_properties(body, properties));
        REQUIRE(f.state.body_pool.set_constant_force(body, {properties.mass * 2, 0}, 0));
        const sqp::tick_t span = sqp::config::TU_PER_SECOND / 3 + 1;
        f.state.interval.target_end = span;
        REQUIRE(sqp::run_pre_solve(f.state));
        REQUIRE(sqp::run_solve(f.state));
        const sqp::real seconds = sqp::math::tu_to_seconds<sqp::real>(span);
        CHECK(f.body(body).linear_velocity.x == Catch::Approx(2 * seconds));
        CHECK(f.body(body).position.x == Catch::Approx(seconds * seconds * (steps + 1) / steps));
        CHECK(f.state.island_pool.get(id)->time == span);
        REQUIRE(f.state.solve.processed_islands.size() == 1);
        CHECK(f.state.solve.processed_islands[0].start == 0);
        CHECK(f.state.solve.processed_islands[0].end == span);
    }
}

TEST_CASE("Contact-connected islands synchronize before ordinary substeps", "[world][solve]") {
    Fixture f;
    f.state.settings.gravity = {};
    f.state.settings.contact.substeps = 4;
    const auto early = f.rigid();
    const auto late = f.rigid({0, 1});
    const auto independent = f.rigid({10, 0});
    const auto early_id = f.island({early});
    const auto late_id = f.island({late}, sqp::config::TU_PER_SECOND / 4);
    const auto independent_id = f.island({independent}, sqp::config::TU_PER_SECOND * 3 / 8);
    sqp::BodyProperties properties;
    REQUIRE(f.state.body_pool.get_properties(early, properties));
    REQUIRE(f.state.body_pool.set_constant_force(early, {properties.mass * 2, 0}, 0));
    sqp::PhysicsMaterial material;
    material.static_friction = 0;
    material.dynamic_friction = 0;
    REQUIRE(f.state.shape_pool.materials().set(0, material));
    f.contact(early, late);
    f.state.interval.target_end = sqp::config::TU_PER_SECOND / 2;
    REQUIRE(sqp::run_pre_solve(f.state));
    REQUIRE(sqp::run_solve(f.state));
    CHECK(f.body(early).position.x == Catch::Approx(0.328125));
    CHECK(f.body(late).position.x == Catch::Approx(0));
    CHECK(f.body(independent).position.x == Catch::Approx(10));
    CHECK(f.state.island_pool.get(early_id)->time == f.state.interval.target_end);
    CHECK(f.state.island_pool.get(late_id)->time == f.state.interval.target_end);
    CHECK(f.state.island_pool.get(independent_id)->time == f.state.interval.target_end);
    CHECK(f.state.solve.processed_islands.size() == 3);
}

TEST_CASE("Static contact endpoints do not synchronize unrelated islands", "[world][solve]") {
    Fixture f;
    f.state.settings.gravity = {};
    f.state.settings.contact.substeps = 4;
    const auto fixed = f.rigid({}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto early = f.rigid({0, 1});
    const auto late = f.rigid({0, 1});
    f.island({fixed}, sqp::config::TU_PER_SECOND * 3 / 8);
    f.island({early});
    f.island({late}, sqp::config::TU_PER_SECOND / 4);
    sqp::BodyProperties properties;
    REQUIRE(f.state.body_pool.get_properties(early, properties));
    REQUIRE(f.state.body_pool.set_constant_force(early, {properties.mass * 2, 0}, 0));
    sqp::PhysicsMaterial material;
    material.static_friction = 0;
    material.dynamic_friction = 0;
    REQUIRE(f.state.shape_pool.materials().set(0, material));
    f.contact(fixed, early);
    f.contact(fixed, late);
    f.state.interval.target_end = sqp::config::TU_PER_SECOND / 2;
    REQUIRE(sqp::run_pre_solve(f.state));
    REQUIRE(sqp::run_solve(f.state));
    CHECK(f.body(early).position.x == Catch::Approx(0.3125));
    CHECK(f.body(late).position.x == Catch::Approx(0));
}

TEST_CASE("TOI microstep remains one integration step", "[world][solve]") {
    Fixture f;
    f.state.settings.gravity = {};
    f.state.settings.contact.substeps = 16;
    const auto body = f.rigid();
    const auto fixed = f.rigid({0, 1}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto id = f.island({body}, 0, false);
    sqp::BodyProperties properties;
    REQUIRE(f.state.body_pool.get_properties(body, properties));
    REQUIRE(f.state.body_pool.set_constant_force(body, {properties.mass * 2, 0}, 0));
    f.contact(body, fixed);
    f.state.interval.microstep = true;
    f.state.interval.microstep_time = sqp::config::TU_PER_SECOND / 2;
    REQUIRE(sqp::run_pre_solve(f.state));
    f.state.pairs.manifolds->counts[0] = 0;
    REQUIRE(sqp::run_solve(f.state));
    CHECK(f.body(body).linear_velocity.x == Catch::Approx(1));
    CHECK(f.body(body).position.x == Catch::Approx(0.5));
    CHECK(f.state.island_pool.get(id)->time == f.state.interval.microstep_time);
}

TEST_CASE("World contact cache carries impulses without positional friction across ordinary and TOI steps",
          "[world][solve][friction][cache]") {
    for (const bool microstep : {false, true}) {
        Fixture f;
        f.state.settings.gravity = {0, -10};
        f.state.settings.contact.substeps = 1;
        const auto ground =
            f.rigid({}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
        const auto body = f.rigid({0, 1}, {}, sqp::BodyConfigs::RegisterType::Dynamic,
                                  sqp::BodyConfigs::OverrideFlags::LockRot);
        f.island({ground});
        f.island({body});
        sqp::PhysicsMaterial material;
        material.static_friction = 4;
        material.dynamic_friction = 2;
        REQUIRE(f.state.shape_pool.materials().set(0, material));
        const auto span = sqp::config::TU_PER_SECOND / 10;
        f.contact(ground, body);
        f.state.interval.target_end = span;
        REQUIRE(sqp::run_pre_solve(f.state));
        REQUIRE(sqp::run_solve(f.state));
        REQUIRE(f.state.solve.contact_cache.size() == 1);
        REQUIRE(f.state.solve.contact_cache[0].count == 1);
        const auto cached_normal = f.state.solve.contact_cache[0].normal_impulses[0];
        CHECK(cached_normal > 0);
        CHECK(f.state.solve.contact_cache[0].step_seconds == Catch::Approx(0.1));

        auto moved = f.body(body);
        moved.position.x += sqp::real(0.01);
        moved.linear_velocity.x = 0;
        moved.linear_velocity.y = -1;
        REQUIRE(f.state.body_pool.set_state(body, moved));
        f.state.pairs.pairs->clear();
        f.state.pairs.manifolds->clear();
        f.state.pairs.original_manifolds.clear();
        f.state.pairs.original_manifolds_valid = false;
        f.contact(ground, body);
        if (microstep) {
            f.state.interval.microstep = true;
            f.state.interval.microstep_time = span * 2;
            f.state.pairs.toi_local_contacts.push_back(
                {{0, 0}, f.state.pairs.manifolds->at(0), true});
        } else {
            f.state.interval.target_end = span * 2;
        }
        REQUIRE(sqp::run_pre_solve(f.state));
        REQUIRE(sqp::run_solve(f.state));
        CHECK(f.body(body).linear_velocity.x == Catch::Approx(0));
        REQUIRE(f.state.solve.contact_cache.size() == 1);
        CHECK(f.state.solve.contact_cache[0].time == span * 2);
    }
}

TEST_CASE("World warm starts a matched contact impulse and includes it in reports",
          "[world][solve][warm][report]") {
    Fixture f;
    f.state.settings.gravity = {0, -10};
    f.state.settings.contact.substeps = 1;
    const auto ground = f.rigid({}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto body = f.rigid({0, 1}, {}, sqp::BodyConfigs::RegisterType::Dynamic,
                              sqp::BodyConfigs::OverrideFlags::LockRot);
    f.island({ground});
    f.island({body});
    sqp::PhysicsMaterial material;
    material.static_friction = 4;
    REQUIRE(f.state.shape_pool.materials().set(0, material));
    REQUIRE(f.state.body_pool.set_event_flags(body, sqp::BodyConfigs::EventFlags::AfterSolve));
    const auto span = sqp::config::TU_PER_SECOND / 10;
    f.contact(ground, body);
    f.state.interval.target_end = span;
    REQUIRE(sqp::run_pre_solve(f.state));
    REQUIRE(sqp::run_solve(f.state));
    REQUIRE(f.state.solve.contact_cache.size() == 1);
    const auto cached_normal = f.state.solve.contact_cache[0].normal_impulses[0];
    REQUIRE(cached_normal > 0);
    const auto seeded_tangent = cached_normal * sqp::real(0.25);
    f.state.solve.contact_cache[0].tangent_impulses[0] = seeded_tangent;
    auto state = f.body(body);
    state.linear_velocity = {0, -1};
    REQUIRE(f.state.body_pool.set_state(body, state));
    f.state.pairs.pairs->clear();
    f.state.pairs.manifolds->clear();
    f.state.pairs.original_manifolds.clear();
    f.state.pairs.original_manifolds_valid = false;
    f.contact(ground, body);
    f.state.interval.target_end = span * 2;
    REQUIRE(sqp::run_pre_solve(f.state));
    REQUIRE(sqp::run_solve(f.state));
    REQUIRE(f.state.post_solve.reports.size() == 1);
    CHECK(f.state.post_solve.reports[0].contact.tangent_impulse >= seeded_tangent);
}

TEST_CASE("World contact cache rejects a changed feature and an interrupted contact",
          "[world][solve][friction][anchor]") {
    for (const bool changed_feature : {false, true}) {
        Fixture f;
        f.state.settings.gravity = {0, -10};
        f.state.settings.contact.substeps = 1;
        const auto ground =
            f.rigid({}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
        const auto body = f.rigid({0, 1}, {}, sqp::BodyConfigs::RegisterType::Dynamic,
                                  sqp::BodyConfigs::OverrideFlags::LockRot);
        f.island({ground});
        f.island({body});
        sqp::PhysicsMaterial material;
        material.static_friction = 4;
        REQUIRE(f.state.shape_pool.materials().set(0, material));
        const auto span = sqp::config::TU_PER_SECOND / 10;
        f.contact(ground, body);
        f.state.interval.target_end = span;
        REQUIRE(sqp::run_pre_solve(f.state));
        REQUIRE(sqp::run_solve(f.state));
        REQUIRE(f.state.solve.contact_cache.size() == 1);

        auto moved = f.body(body);
        moved.position.x += sqp::real(0.01);
        moved.linear_velocity.x = 0;
        REQUIRE(f.state.body_pool.set_state(body, moved));
        f.state.pairs.pairs->clear();
        f.state.pairs.manifolds->clear();
        f.state.pairs.original_manifolds.clear();
        f.state.pairs.original_manifolds_valid = false;
        if (changed_feature) {
            f.contact(ground, body);
            f.state.pairs.manifolds->keys[0].feature_a.index = 1;
        }
        f.state.interval.target_end = span * 2;
        REQUIRE(sqp::run_pre_solve(f.state));
        REQUIRE(sqp::run_solve(f.state));
        if (!changed_feature) {
            CHECK(f.state.solve.contact_cache.empty());
            f.state.pairs.pairs->clear();
            f.state.pairs.manifolds->clear();
            f.state.pairs.original_manifolds.clear();
            f.state.pairs.original_manifolds_valid = false;
            f.contact(ground, body);
            f.state.interval.target_end = span * 3;
            REQUIRE(sqp::run_pre_solve(f.state));
            REQUIRE(sqp::run_solve(f.state));
        }
        CHECK(f.body(body).linear_velocity.x == Catch::Approx(0).margin(1e-5));
    }
}

TEST_CASE("World contact cache does not follow a recycled shape slot",
          "[world][solve][friction][anchor]") {
    Fixture f;
    f.state.settings.gravity = {0, -10};
    f.state.settings.contact.substeps = 1;
    const auto ground = f.rigid({}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto body = f.rigid({0, 1}, {}, sqp::BodyConfigs::RegisterType::Dynamic,
                              sqp::BodyConfigs::OverrideFlags::LockRot);
    f.island({ground});
    f.island({body});
    sqp::PhysicsMaterial material;
    material.static_friction = 4;
    REQUIRE(f.state.shape_pool.materials().set(0, material));
    const auto span = sqp::config::TU_PER_SECOND / 10;
    f.contact(ground, body);
    f.state.interval.target_end = span;
    REQUIRE(sqp::run_pre_solve(f.state));
    REQUIRE(sqp::run_solve(f.state));
    REQUIRE(f.state.solve.contact_cache.size() == 1);
    const auto old_shape = f.state.solve.contact_cache[0].shape_b;
    REQUIRE(f.state.body_pool.detach_shape(body, old_shape, f.state.shape_pool));
    REQUIRE(f.state.shape_pool.remove(old_shape));
    const auto new_shape = f.state.shape_pool.add_circle({sqp::real(1)});
    REQUIRE(f.state.body_pool.attach_shape(body, new_shape, f.state.shape_pool));
    CHECK(new_shape != old_shape);

    auto moved = f.body(body);
    moved.position.x += sqp::real(0.01);
    moved.linear_velocity.x = 0;
    REQUIRE(f.state.body_pool.set_state(body, moved));
    f.state.pairs.pairs->clear();
    f.state.pairs.manifolds->clear();
    f.state.pairs.original_manifolds.clear();
    f.state.pairs.original_manifolds_valid = false;
    f.contact(ground, body);
    f.state.interval.target_end = span * 2;
    REQUIRE(sqp::run_pre_solve(f.state));
    REQUIRE(sqp::run_solve(f.state));
    CHECK(f.body(body).linear_velocity.x == Catch::Approx(0).margin(1e-5));
    REQUIRE(f.state.solve.contact_cache.size() == 1);
    CHECK(f.state.solve.contact_cache[0].shape_b == new_shape);
}

TEST_CASE("World TOI rebounds a circle from a wall near x=1024 in either pair order",
          "[world][solve][toi]") {
    for (const bool reverse : {false, true}) {
        Fixture f;
        f.state.settings.gravity = {};
        f.state.interval.target_end = sqp::config::TU_PER_SECOND / 10;
        sqp::RigidBodyDef wall_def;
        wall_def.position = {sqp::real(1024.5), 0};
        wall_def.configs.register_type = sqp::BodyConfigs::RegisterType::ImmutableStatic;
        const auto wall = f.state.body_pool.add_rigid(wall_def);
        sqp::PolygonShapeDef wall_shape_def;
        wall_shape_def.vertices = {{sqp::real(-0.5), -10}, {sqp::real(0.5), -10},
                                   {sqp::real(0.5), 10}, {sqp::real(-0.5), 10}};
        const auto wall_shape = f.state.shape_pool.add_polygon(wall_shape_def);
        REQUIRE(f.state.body_pool.attach_shape(wall, wall_shape, f.state.shape_pool));
        const auto moving = f.rigid({sqp::real(1022.5), 0}, {10, 0});
        f.island({moving}, 0, false);
        std::vector<sqp::ShapeHandle> shapes;
        REQUIRE(f.state.body_pool.get_attached_shapes(moving, shapes));
        REQUIRE(shapes.size() == 1);
        const auto circle_shape = shapes[0];
        if (reverse) {
            f.state.broad_phase.bodies = {moving, wall};
            f.state.broad_phase.shapes = {circle_shape, wall_shape};
        } else {
            f.state.broad_phase.bodies = {wall, moving};
            f.state.broad_phase.shapes = {wall_shape, circle_shape};
        }
        f.state.broad_phase.candidates.push_back({0, 1});
        sqp::PhysicsMaterial material;
        material.base_restitution = 1;
        material.restitution_flags = sqp::RESTITUTION_DISABLE_IMPACT_SCALING;
        REQUIRE(f.state.shape_pool.materials().set(0, material));
        REQUIRE(sqp::run_narrow_phase(f.state));
        REQUIRE(f.state.narrow_phase.toi_pairs.size() == 1);
        REQUIRE(f.state.narrow_phase.toi_impacts.size() == 1);
        bool selected = false;
        REQUIRE(sqp::run_process_toi(f.state, f.state.interval.target_end, selected));
        REQUIRE(selected);
        REQUIRE(f.state.pairs.toi_local_contacts.size() == 1);
        CHECK(f.state.pairs.toi_local_contacts[0].valid);
        REQUIRE(sqp::run_pre_solve(f.state));
        REQUIRE(sqp::run_solve(f.state));
        CHECK(f.body(moving).linear_velocity.x < 0);
    }
}

TEST_CASE("Consecutive TOI impulses preserve the position reached before the second impact",
          "[world][solve][toi]") {
    Fixture f;
    f.state.settings.gravity = {};
    const auto moving = f.rigid({}, {0, 4});
    const auto upper = f.rigid({0, 2}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto lower = f.rigid({0, -2}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto island = f.island({moving}, 0, false);
    f.contact(moving, upper);
    f.contact(lower, moving);
    f.state.pairs.manifolds->counts[1] = 0;
    sqp::PhysicsMaterial material;
    material.base_restitution = 1;
    material.restitution_flags = sqp::RESTITUTION_DISABLE_IMPACT_SCALING;
    material.static_friction = 0;
    material.dynamic_friction = 0;
    REQUIRE(f.state.shape_pool.materials().set(0, material));

    f.state.interval.microstep = true;
    f.state.interval.microstep_time = sqp::config::TU_PER_SECOND / 4;
    REQUIRE(sqp::run_pre_solve(f.state));
    REQUIRE(sqp::run_solve(f.state));
    const auto first = f.body(moving);
    REQUIRE(first.linear_velocity.y < 0);
    CHECK(first.position.y == Catch::Approx(1));

    f.state.pairs.manifolds->counts[0] = 0;
    f.state.pairs.manifolds->counts[1] = 1;
    f.state.interval.microstep_time = sqp::config::TU_PER_SECOND / 2;
    REQUIRE(sqp::run_pre_solve(f.state));
    REQUIRE(sqp::run_solve(f.state));
    const auto second = f.body(moving);
    const auto seconds = sqp::math::tu_to_seconds<sqp::real>(sqp::config::TU_PER_SECOND / 4);
    CHECK(second.position.y == Catch::Approx(first.position.y + first.linear_velocity.y * seconds));
    CHECK(second.linear_velocity.y > 0);
    CHECK(f.state.island_pool.get(island)->time == f.state.interval.microstep_time);
}

TEST_CASE("Solve resolves a PreSolve contact and aligns impulse reports", "[world][solve]") {
    Fixture f;
    f.state.settings.gravity = {};
    const auto ground = f.rigid({}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto moving = f.rigid({0, 1}, {0, -2});
    const auto ground_island = f.island({ground});
    const auto moving_island = f.island({moving});
    f.contact(ground, moving);
    f.contact(ground, moving);
    sqp::PhysicsMaterial material;
    material.base_restitution = 1;
    material.restitution_flags = sqp::RESTITUTION_DISABLE_IMPACT_SCALING;
    REQUIRE(f.state.shape_pool.materials().set(0, material));
    REQUIRE(sqp::run_pre_solve(f.state));
    f.state.pairs.manifolds->counts[0] = 0;
    REQUIRE(f.state.body_pool.set_event_flags(moving, sqp::BodyConfigs::EventFlags::AfterSolve));
    f.state.interval.target_end = sqp::config::TU_PER_SECOND / 64;

    REQUIRE(sqp::run_solve(f.state));
    CHECK(f.state.island_pool.get(ground_island)->time == f.state.interval.target_end);
    CHECK(f.state.island_pool.get(moving_island)->time == f.state.interval.target_end);
    CHECK(f.body(ground).position.y == 0);
    CHECK(f.body(moving).linear_velocity.y > -2);
    REQUIRE(f.state.post_solve.reports.size() == 1);
    CHECK(f.state.post_solve.reports[0].pair_index == 1);
    CHECK(f.state.post_solve.reports[0].contact.normal_impulse > 0);
    CHECK(f.state.post_solve.reports[0].contact.points[0].normal_impulse > 0);
    CHECK(f.state.post_solve.reports[0].contact.points[1].normal_impulse == 0);
    CHECK(f.state.solve.constraints.empty());
    CHECK(f.state.pairs.combined_restitutions.empty());
}

TEST_CASE("Solve microstep advances pair islands from their own times", "[world][solve]") {
    Fixture f;
    f.state.settings.gravity = {};
    const auto a = f.rigid({}, {1, 0});
    const auto b = f.rigid({0, 1}, {2, 0});
    const auto early = f.island({a}, 0, false);
    const auto late = f.island({b}, sqp::config::TU_PER_SECOND / 4, false);
    f.contact(a, b);
    REQUIRE(sqp::run_pre_solve(f.state));
    f.state.pairs.manifolds->counts[0] = 0;
    f.state.interval.microstep = true;
    f.state.interval.microstep_time = sqp::config::TU_PER_SECOND / 2;

    REQUIRE(sqp::run_solve(f.state));
    CHECK(f.body(a).position.x == Catch::Approx(0.5));
    CHECK(f.body(b).position.x == Catch::Approx(0.5));
    CHECK(f.state.island_pool.get(early)->time == f.state.interval.microstep_time);
    CHECK(f.state.island_pool.get(late)->time == f.state.interval.microstep_time);
    CHECK(f.state.post_solve.reports.empty());
}

TEST_CASE("Solve omits reports for contacts already at the target time", "[world][solve]") {
    Fixture f;
    f.state.settings.gravity = {};
    const auto a = f.rigid();
    const auto b = f.rigid({0, 1});
    const auto target = sqp::config::TU_PER_SECOND / 4;
    f.island({a, b}, target);
    f.contact(a, b);
    REQUIRE(f.state.body_pool.set_event_flags(b, sqp::BodyConfigs::EventFlags::AfterSolve));
    REQUIRE(sqp::run_pre_solve(f.state));
    f.state.interval.target_end = target;

    REQUIRE(sqp::run_solve(f.state));
    CHECK(f.state.post_solve.reports.empty());
    CHECK(f.state.solve.constraints.empty());
}

TEST_CASE("Solve uses the shared remaining time for a mixed-time contact", "[world][solve]") {
    Fixture f;
    f.state.settings.gravity = {};
    const auto ground = f.rigid({}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto moving = f.rigid({0, 1}, {0, -2});
    const auto ground_island = f.island({ground}, 0, false);
    const auto moving_island = f.island({moving}, sqp::config::TU_PER_SECOND / 4, false);
    f.contact(ground, moving);
    REQUIRE(sqp::run_pre_solve(f.state));
    REQUIRE(f.state.body_pool.set_event_flags(moving, sqp::BodyConfigs::EventFlags::AfterSolve));
    f.state.interval.microstep = true;
    f.state.interval.microstep_time = sqp::config::TU_PER_SECOND / 2;

    REQUIRE(sqp::run_solve(f.state));
    REQUIRE(f.state.post_solve.reports.size() == 1);
    CHECK(f.state.post_solve.reports[0].pair_index == 0);
    CHECK(f.state.post_solve.reports[0].contact.normal_impulse > 0);
    CHECK(f.state.island_pool.get(ground_island)->time == f.state.interval.microstep_time);
    CHECK(f.state.island_pool.get(moving_island)->time == f.state.interval.microstep_time);
}

TEST_CASE("Solve keeps kinematic and static body behavior and locked axes", "[world][solve]") {
    Fixture f;
    f.state.settings.gravity = {3, -10};
    const auto kinematic = f.rigid({}, {2, 3}, sqp::BodyConfigs::RegisterType::Kinematic);
    const auto fixed = f.rigid({1, 2}, {3, 4}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto locked = f.rigid({}, {}, sqp::BodyConfigs::RegisterType::Dynamic,
                                sqp::BodyConfigs::OverrideFlags::LockLinearX |
                                    sqp::BodyConfigs::OverrideFlags::LockRot);
    const auto no_gravity_y = f.rigid({}, {}, sqp::BodyConfigs::RegisterType::Dynamic,
                                      sqp::BodyConfigs::OverrideFlags::NoGravityY);
    f.island({kinematic, fixed, locked, no_gravity_y});
    REQUIRE(f.state.body_pool.set_constant_force(locked, {10, 0}, 5));
    f.state.interval.target_end = sqp::config::TU_PER_SECOND / 2;
    REQUIRE(sqp::run_pre_solve(f.state));
    REQUIRE(sqp::run_solve(f.state));
    CHECK(f.body(kinematic).position == sqp::Vec2{1, sqp::real(1.5)});
    CHECK(f.body(kinematic).linear_velocity == sqp::Vec2{2, 3});
    CHECK(f.body(fixed).position == sqp::Vec2{1, 2});
    CHECK(f.body(locked).position.x == 0);
    CHECK(f.body(locked).linear_velocity.x == 0);
    CHECK(f.body(locked).angular_velocity == 0);
    CHECK(f.body(locked).linear_velocity.y == Catch::Approx(-5));
    CHECK(f.body(no_gravity_y).linear_velocity.y == 0);
    CHECK(f.body(no_gravity_y).linear_velocity.x == Catch::Approx(1.5));
}

TEST_CASE("Solve rejects invalid input before advancing bodies", "[world][solve]") {
    Fixture f;
    const auto body = f.rigid();
    const auto id = f.island({body});
    f.state.interval.target_end = sqp::config::TU_PER_SECOND / 2;
    REQUIRE(sqp::run_pre_solve(f.state));
    f.state.interval.target_end = 0;
    REQUIRE(f.state.island_pool.set_time(id, 1));
    const auto before = f.body(body);
    CHECK_FALSE(sqp::run_solve(f.state));
    CHECK(f.body(body) == before);
    CHECK(f.state.island_pool.get(id)->time == 1);
}

TEST_CASE("Solve rejects missing PreSolve output without changing a body", "[world][solve]") {
    Fixture f;
    const auto a = f.rigid();
    const auto b = f.rigid({0, 1});
    const auto id = f.island({a, b});
    f.contact(a, b);
    REQUIRE(sqp::run_pre_solve(f.state));
    f.state.interval.target_end = sqp::config::TU_PER_SECOND / 4;
    f.state.pairs.combined_dynamic_frictions.clear();
    const auto before = f.body(b);

    CHECK_FALSE(sqp::run_solve(f.state));
    CHECK(f.body(b) == before);
    CHECK(f.state.island_pool.get(id)->time == 0);
}

TEST_CASE("Solve collects mixed body kinds and repeated island references once", "[world][solve]") {
    Fixture f;
    f.state.settings.gravity = {};
    const auto first = f.rigid({}, {1, 0});
    const auto second = f.rigid({}, {2, 0});
    const auto circle_shape = f.state.shape_pool.add_circle({1});
    const auto circle = f.state.body_pool.add_lite_circle({}, circle_shape, f.state.shape_pool);
    sqp::PolygonShapeDef box_def;
    box_def.vertices = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
    const auto box_shape = f.state.shape_pool.add_polygon(box_def);
    const auto box = f.state.body_pool.add_aabb_box({}, box_shape, f.state.shape_pool);
    const auto first_island = f.island({first, circle, box});
    const auto second_island = f.island({second});
    f.state.interval.islands->push_back(first_island);
    f.contact(first, second);
    f.contact(first, second);
    REQUIRE(sqp::run_pre_solve(f.state));
    f.state.pairs.manifolds->counts[0] = 0;
    f.state.pairs.manifolds->counts[1] = 0;
    f.state.interval.target_end = sqp::config::TU_PER_SECOND / 4;

    REQUIRE(sqp::run_solve(f.state));
    CHECK(f.state.solve.islands.empty());
    CHECK(f.state.solve.bodies.empty());
    CHECK(f.state.solve.body_start_times.empty());
    CHECK(f.state.post_solve.reports.empty());
    CHECK(f.state.island_pool.get(first_island)->time == f.state.interval.target_end);
    CHECK(f.state.island_pool.get(second_island)->time == f.state.interval.target_end);
    CHECK(f.body(first).position.x == Catch::Approx(0.25));
    CHECK(f.body(second).position.x == Catch::Approx(0.5));
}

TEST_CASE("Solve stamps survive consecutive calls, wraparound and reused island slots",
          "[world][solve]") {
    Fixture f;
    f.state.settings.gravity = {};
    const auto body = f.rigid({}, {1, 0});
    const auto old_id = f.island({body});
    REQUIRE(sqp::run_pre_solve(f.state));
    f.state.interval.target_end = sqp::config::TU_PER_SECOND / 4;
    REQUIRE(sqp::run_solve(f.state));

    f.state.solve.selection_epoch = std::numeric_limits<sqp::u32>::max();
    f.state.interval.target_end = sqp::config::TU_PER_SECOND / 2;
    REQUIRE(sqp::run_solve(f.state));
    CHECK(f.state.solve.selection_epoch == 1);
    CHECK(f.state.solve.bodies.empty());
    CHECK(f.body(body).position.x == Catch::Approx(0.5));

    REQUIRE(f.state.island_pool.remove_body(old_id, body, f.state.body_pool));
    REQUIRE(f.state.island_pool.destroy(old_id));
    const auto new_id = f.state.island_pool.create();
    REQUIRE(new_id.index == old_id.index);
    REQUIRE(new_id.generation != old_id.generation);
    REQUIRE(f.state.island_pool.set_time(new_id, sqp::config::TU_PER_SECOND / 2));
    REQUIRE(f.state.island_pool.add_body(new_id, body, f.state.body_pool));
    f.state.interval.islands->assign(1, new_id);
    f.state.interval.target_end = sqp::config::TU_PER_SECOND * 3 / 4;
    REQUIRE(sqp::run_solve(f.state));
    CHECK(f.state.solve.islands.empty());
    CHECK(f.state.solve.bodies.empty());
    CHECK(f.body(body).position.x == Catch::Approx(0.75));

    const auto old_dense = f.state.body_pool.dense_index(body);
    REQUIRE(old_dense);
    REQUIRE(f.state.body_pool.remove(body, f.state.island_pool));
    const auto replacement = f.rigid({}, {2, 0});
    CHECK(f.state.body_pool.dense_index(replacement) == old_dense);
    REQUIRE(f.state.island_pool.add_body(new_id, replacement, f.state.body_pool));
    f.state.interval.target_end = sqp::config::TU_PER_SECOND;
    REQUIRE(sqp::run_solve(f.state));
    CHECK(f.state.solve.bodies.empty());
    CHECK(f.body(replacement).position.x == Catch::Approx(0.5));
}

#ifdef SQP_DEBUG_VALIDATION
TEST_CASE("Solve rejects duplicate body membership before integration", "[world][solve]") {
    Fixture f;
    const auto body = f.rigid({}, {1, 0});
    const auto id = f.island({body});
    auto* island = const_cast<sqp::Island*>(f.state.island_pool.get(id));
    REQUIRE(island);
    std::get<sqp::Island::SmallIsland>(island->bodies).data[1] = body;
    REQUIRE(sqp::run_pre_solve(f.state));
    f.state.interval.target_end = sqp::config::TU_PER_SECOND / 4;
    const auto before = f.body(body);

    CHECK_FALSE(sqp::run_solve(f.state));
    CHECK(f.body(body) == before);
    CHECK(f.state.island_pool.get(id)->time == 0);
}
#endif

TEST_CASE("Solve collects many active islands and bodies", "[world][solve]") {
    Fixture f;
    f.state.settings.gravity = {};
    constexpr std::size_t count = 512;
    for (std::size_t i = 0; i < count; ++i) {
        const auto body = f.state.body_pool.add_rigid();
        REQUIRE(f.state.body_pool.is_valid(body));
        f.island({body});
    }
    REQUIRE(sqp::run_pre_solve(f.state));
    f.state.interval.target_end = sqp::config::TU_PER_SECOND / 4;

    REQUIRE(sqp::run_solve(f.state));
    CHECK(f.state.solve.islands.empty());
    CHECK(f.state.solve.bodies.empty());
    CHECK(f.state.solve.body_start_times.empty());
    CHECK(f.state.post_solve.reports.empty());
    CHECK(f.state.island_pool.get(f.state.interval.islands->back())->time ==
          f.state.interval.target_end);
}
