#include "world/world_stages.h"

#include <catch_amalgamated.hpp>

#include <cmath>
#include <limits>

namespace {

struct Fixture {
    sqp::WorldState state;
    sqp::BodyHandle a, b;
    sqp::ShapeHandle shape_a, shape_b;

    Fixture() {
        a = state.body_pool.add_rigid({{}, 0, {}, 0});
        b = state.body_pool.add_rigid({{0, 1}, 0, {0, -2}, 0});
        shape_a = state.shape_pool.add_circle({sqp::real(0.5)});
        shape_b = state.shape_pool.add_circle({sqp::real(0.5)});
        REQUIRE(state.body_pool.attach_shape(a, shape_a, state.shape_pool));
        REQUIRE(state.body_pool.attach_shape(b, shape_b, state.shape_pool));
    }

    void append() {
        const auto versions_a = state.body_pool.motion_version_column(a.kind());
        const auto versions_b = state.body_pool.motion_version_column(b.kind());
        state.pairs.pairs->push_back({a, b, shape_a, shape_b,
                                      versions_a[*state.body_pool.dense_index(a)],
                                      versions_b[*state.body_pool.dense_index(b)]});
        sqp::ContactManifold manifold;
        manifold.count = 1;
        manifold.normal = {0, 1};
        manifold.points[0].point = {0, sqp::real(0.5)};
        manifold.points[0].depth = sqp::real(0.01);
        REQUIRE(state.pairs.manifolds->append(manifold));
    }
};

struct CallbackData {
    int calls = 0;
    enum class Action {
        Adjust,
        Disable,
        InvalidNumber,
        InvalidPair,
        InvalidNormal,
        SecondInvalid,
        Rows,
        MaterialOnly
    } action = Action::Adjust;
};

void callback(sqp::PhysicsEvent& event, void* context) {
    auto& data = *static_cast<CallbackData*>(context);
    ++data.calls;
    auto& contact = std::get<sqp::PreSolveEvent>(event.event);
    switch (data.action) {
    case CallbackData::Action::Adjust:
        contact.manifold.points[0].point.x = sqp::real(0.25);
        contact.combined_restitution = sqp::real(0.6);
        contact.combined_dynamic_friction = sqp::real(0.2);
        contact.combined_static_friction = sqp::real(0.4);
        contact.combined_damping_ratio = sqp::real(3);
        contact.combined_contact_frequency = sqp::real(12);
        contact.combined_max_penetration_correction_speed = sqp::real(2);
        break;
    case CallbackData::Action::Disable:
        contact.manifold.count = 0;
        break;
    case CallbackData::Action::InvalidNumber:
        contact.combined_restitution = std::numeric_limits<sqp::real>::quiet_NaN();
        break;
    case CallbackData::Action::InvalidPair:
        contact.contact_pair.body_a = contact.contact_pair.body_b;
        break;
    case CallbackData::Action::InvalidNormal:
        contact.manifold.normal = {};
        break;
    case CallbackData::Action::SecondInvalid:
        if (data.calls == 1)
            contact.manifold.points[0].point.x = sqp::real(0.25);
        else
            contact.combined_dynamic_friction = -1;
        break;
    case CallbackData::Action::Rows:
        contact.manifold.points[0].point.x = sqp::real(data.calls) / sqp::real(4);
        contact.combined_restitution = sqp::real(data.calls) / sqp::real(10);
        contact.combined_dynamic_friction = sqp::real(data.calls) / sqp::real(10);
        contact.combined_static_friction = sqp::real(data.calls) / sqp::real(5);
        contact.combined_damping_ratio = sqp::real(data.calls);
        contact.combined_contact_frequency = sqp::real(data.calls * 10);
        contact.combined_max_penetration_correction_speed = sqp::real(data.calls * 2);
        if (data.calls == 1)
            contact.manifold.count = 0;
        break;
    case CallbackData::Action::MaterialOnly:
        contact.combined_restitution = sqp::real(0.5);
        break;
    }
}

} // namespace

TEST_CASE("PreSolve preserves TOI local geometry for material edits and remaps geometry edits",
          "[world][pre_solve][toi]") {
    for (auto action : {CallbackData::Action::MaterialOnly, CallbackData::Action::Adjust}) {
        Fixture f;
        f.append();
        const sqp::Vec2 origin{1024, 0};
        auto local = f.state.pairs.manifolds->at(0);
        local.points[0].point = {sqp::real(0.123456), sqp::real(0.5)};
        f.state.pairs.toi_local_contacts.push_back({origin, local, true});
        CallbackData data;
        data.action = action;
        f.state.callbacks = {callback, nullptr, &data, nullptr};
        REQUIRE(f.state.body_pool.set_event_flags(f.a, sqp::BodyConfigs::EventFlags::PreSolve));
        REQUIRE(sqp::run_pre_solve(f.state));
        REQUIRE(data.calls == 1);
        const auto& retained = f.state.pairs.toi_local_contacts[0].manifold;
        if (action == CallbackData::Action::MaterialOnly)
            CHECK(retained.points[0].point.x == sqp::real(0.123456));
        else
            CHECK(retained.points[0].point.x == sqp::real(0.25) - origin.x);
    }
}

TEST_CASE("PreSolve accepts empty contacts and refreshes output columns", "[world][pre_solve]") {
    Fixture f;
    f.state.pairs.combined_restitutions.push_back(1);
    REQUIRE(sqp::run_pre_solve(f.state));
    CHECK(f.state.pairs.combined_restitutions.empty());
    CHECK(f.state.pairs.combined_static_frictions.empty());
}

TEST_CASE("PreSolve computes ordinary and microstep restitution", "[world][pre_solve]") {
    for (bool microstep : {false, true}) {
        Fixture f;
        f.append();
        f.state.interval.microstep = microstep;
        REQUIRE(sqp::run_pre_solve(f.state));
        REQUIRE(f.state.pairs.combined_restitutions.size() == 1);
        CHECK(f.state.pairs.combined_restitutions[0] >= 0);
        CHECK(f.state.pairs.combined_restitutions[0] <= 1);
        CHECK(f.state.pairs.combined_static_frictions[0] >=
              f.state.pairs.combined_dynamic_frictions[0]);
    }
}

TEST_CASE("PreSolve dynamic restitution uses the lesser shape impact scale", "[world][pre_solve]") {
    Fixture f;
    f.append();
    sqp::PhysicsMaterial material;
    material.base_restitution = 1;
    material.restitution_impact_threshold = 1;
    material.restitution_impact_exponent = 1;
    REQUIRE(f.state.shape_pool.materials().set(0, material));
    REQUIRE(f.state.shape_pool.set_restitution_impact_scale(f.shape_a, 2));
    REQUIRE(f.state.shape_pool.set_restitution_impact_scale(f.shape_b, 4));
    REQUIRE(sqp::run_pre_solve(f.state));
    const auto first = f.state.pairs.combined_restitutions[0];
    REQUIRE(f.state.shape_pool.set_restitution_impact_scale(f.shape_a, 4));
    REQUIRE(sqp::run_pre_solve(f.state));
    CHECK(first > f.state.pairs.combined_restitutions[0]);
}

TEST_CASE("PreSolve has finite restitution without a closing impact", "[world][pre_solve]") {
    for (const sqp::real exponent :
         {sqp::real(0), sqp::real(0.25), sqp::real(0.5), sqp::real(1), sqp::real(2)}) {
        for (const sqp::real velocity : {sqp::real(0), sqp::real(2)}) {
            Fixture f;
            f.append();
            sqp::PhysicsMaterial material;
            material.base_restitution = sqp::real(0.8);
            material.restitution_impact_exponent = exponent;
            REQUIRE(f.state.shape_pool.materials().set(0, material));
            sqp::BodyState body;
            REQUIRE(f.state.body_pool.get_state(f.b, body));
            body.linear_velocity.y = velocity;
            REQUIRE(f.state.body_pool.set_state(f.b, body));

            REQUIRE(sqp::run_pre_solve(f.state));
            REQUIRE(f.state.pairs.combined_restitutions.size() == 1);
            CHECK(std::isfinite(f.state.pairs.combined_restitutions[0]));
            CHECK(f.state.pairs.combined_restitutions[0] == sqp::real(0));
        }
    }
}

TEST_CASE("PreSolve combines shape materials by priority", "[world][pre_solve]") {
    Fixture f;
    f.append();
    sqp::PhysicsMaterial first, second;
    first.dynamic_friction = sqp::real(0.6);
    first.static_friction = sqp::real(0.8);
    second.dynamic_friction = sqp::real(0.2);
    second.static_friction = sqp::real(0.4);
    const auto first_index = f.state.shape_pool.materials().add(first);
    const auto second_index = f.state.shape_pool.materials().add(second);
    REQUIRE(f.state.shape_pool.set_material_index(f.shape_a, first_index));
    REQUIRE(f.state.shape_pool.set_material_index(f.shape_b, second_index));
    REQUIRE(sqp::run_pre_solve(f.state));
    CHECK(f.state.pairs.combined_dynamic_frictions[0] == sqp::real(0.2));
    first.friction_priority = 1;
    REQUIRE(f.state.shape_pool.materials().set(first_index, first));
    REQUIRE(sqp::run_pre_solve(f.state));
    CHECK(f.state.pairs.combined_dynamic_frictions[0] == sqp::real(0.6));
    CHECK(f.state.pairs.combined_static_frictions[0] == sqp::real(0.8));
}

TEST_CASE("PreSolve callback requires a body flag and writes all overrides", "[world][pre_solve]") {
    Fixture f;
    f.append();
    CallbackData data;
    f.state.callbacks = {callback, nullptr, &data, nullptr};
    REQUIRE(sqp::run_pre_solve(f.state));
    CHECK(data.calls == 0);
    REQUIRE(f.state.body_pool.set_event_flags(f.b, sqp::BodyConfigs::EventFlags::PreSolve));
    REQUIRE(sqp::run_pre_solve(f.state));
    CHECK(data.calls == 1);
    CHECK(f.state.pairs.manifolds->at(0).points[0].point.x == sqp::real(0.25));
    CHECK(f.state.pairs.combined_restitutions[0] == sqp::real(0.6));
    CHECK(f.state.pairs.combined_dynamic_frictions[0] == sqp::real(0.2));
    CHECK(f.state.pairs.combined_static_frictions[0] == sqp::real(0.4));
    CHECK(f.state.pairs.combined_damping_ratios[0] == sqp::real(3));
    CHECK(f.state.pairs.combined_contact_frequencies[0] == sqp::real(12));
    CHECK(f.state.pairs.combined_max_penetration_correction_speeds[0] == sqp::real(2));
}

TEST_CASE("PreSolve disables contact using zero manifold points", "[world][pre_solve]") {
    Fixture f;
    f.append();
    CallbackData data;
    data.action = CallbackData::Action::Disable;
    f.state.callbacks = {callback, nullptr, &data, nullptr};
    REQUIRE(f.state.body_pool.set_event_flags(f.a, sqp::BodyConfigs::EventFlags::PreSolve));
    REQUIRE(sqp::run_pre_solve(f.state));
    CHECK(f.state.pairs.pairs->size() == 1);
    CHECK(f.state.pairs.manifolds->counts[0] == 0);
    CHECK(f.state.pairs.manifolds->keys.size() == 2);
    REQUIRE(sqp::run_pre_solve(f.state));
    CHECK(f.state.pairs.manifolds->counts[0] == 0);
}

TEST_CASE("PreSolve keeps callback output aligned across contact rows", "[world][pre_solve]") {
    Fixture f;
    f.append();
    f.append();
    CallbackData data;
    data.action = CallbackData::Action::Rows;
    f.state.callbacks = {callback, nullptr, &data, nullptr};
    REQUIRE(f.state.body_pool.set_event_flags(f.a, sqp::BodyConfigs::EventFlags::PreSolve));
    REQUIRE(f.state.body_pool.set_event_flags(f.b, sqp::BodyConfigs::EventFlags::PreSolve));

    REQUIRE(sqp::run_pre_solve(f.state));
    CHECK(data.calls == 2);
    REQUIRE(f.state.pairs.pairs->size() == 2);
    REQUIRE(f.state.pairs.manifolds->size() == 2);
    REQUIRE(f.state.pairs.combined_restitutions.size() == 2);
    REQUIRE(f.state.pairs.combined_dynamic_frictions.size() == 2);
    REQUIRE(f.state.pairs.combined_static_frictions.size() == 2);
    REQUIRE(f.state.pairs.combined_damping_ratios.size() == 2);
    REQUIRE(f.state.pairs.combined_contact_frequencies.size() == 2);
    REQUIRE(f.state.pairs.combined_max_penetration_correction_speeds.size() == 2);
    CHECK(f.state.pairs.manifolds->counts[0] == 0);
    CHECK(f.state.pairs.manifolds->counts[1] == 1);
    for (std::size_t i = 0; i < 2; ++i) {
        const sqp::real row = sqp::real(i + 1);
        CHECK(f.state.pairs.manifolds->at(i).points[0].point.x == row / sqp::real(4));
        CHECK(f.state.pairs.combined_restitutions[i] == row / sqp::real(10));
        CHECK(f.state.pairs.combined_dynamic_frictions[i] == row / sqp::real(10));
        CHECK(f.state.pairs.combined_static_frictions[i] == row / sqp::real(5));
        CHECK(f.state.pairs.combined_damping_ratios[i] == row);
        CHECK(f.state.pairs.combined_contact_frequencies[i] == row * sqp::real(10));
        CHECK(f.state.pairs.combined_max_penetration_correction_speeds[i] == row * sqp::real(2));
    }
}

TEST_CASE("PreSolve rejects invalid layout, shape ownership and material", "[world][pre_solve]") {
    Fixture f;
    f.append();
    f.state.pairs.manifolds->depths.pop_back();
    CHECK_FALSE(sqp::run_pre_solve(f.state));
    f.state.pairs.manifolds->depths.push_back(0);
    f.state.pairs.pairs->at(0).shape_a = f.shape_b;
    CHECK_FALSE(sqp::run_pre_solve(f.state));
    f.state.pairs.pairs->at(0).shape_a = f.shape_a;
    sqp::PhysicsMaterial invalid;
    invalid.dynamic_friction = std::numeric_limits<sqp::real>::infinity();
    REQUIRE(f.state.shape_pool.materials().set(0, invalid));
    CHECK_FALSE(sqp::run_pre_solve(f.state));
}

TEST_CASE("PreSolve rejects invalid callbacks without partial writes", "[world][pre_solve]") {
    for (const auto action :
         {CallbackData::Action::InvalidNumber, CallbackData::Action::InvalidPair,
          CallbackData::Action::InvalidNormal}) {
        Fixture f;
        f.append();
        f.append();
        CallbackData data;
        data.action = action;
        f.state.callbacks = {callback, nullptr, &data, nullptr};
        REQUIRE(f.state.body_pool.set_event_flags(f.a, sqp::BodyConfigs::EventFlags::PreSolve));
        const auto old_manifold = f.state.pairs.manifolds->at(0);
        CHECK_FALSE(sqp::run_pre_solve(f.state));
        CHECK(f.state.pairs.combined_restitutions.empty());
        CHECK(f.state.pairs.manifolds->at(0).points[0].point == old_manifold.points[0].point);
    }
}

TEST_CASE("PreSolve leaves prior outputs unchanged after a later callback fails",
          "[world][pre_solve]") {
    Fixture f;
    f.append();
    f.append();
    REQUIRE(sqp::run_pre_solve(f.state));
    const auto old_manifold = f.state.pairs.manifolds->at(0);
    const auto old_restitution = f.state.pairs.combined_restitutions;
    const auto old_dynamic_friction = f.state.pairs.combined_dynamic_frictions;
    const auto old_static_friction = f.state.pairs.combined_static_frictions;
    const auto old_damping = f.state.pairs.combined_damping_ratios;
    const auto old_frequency = f.state.pairs.combined_contact_frequencies;
    const auto old_correction_speed = f.state.pairs.combined_max_penetration_correction_speeds;
    CallbackData data;
    data.action = CallbackData::Action::SecondInvalid;
    f.state.callbacks = {callback, nullptr, &data, nullptr};
    REQUIRE(f.state.body_pool.set_event_flags(f.a, sqp::BodyConfigs::EventFlags::PreSolve));
    CHECK_FALSE(sqp::run_pre_solve(f.state));
    CHECK(data.calls == 2);
    CHECK(f.state.pairs.manifolds->at(0).points[0].point == old_manifold.points[0].point);
    CHECK(f.state.pairs.combined_restitutions == old_restitution);
    CHECK(f.state.pairs.combined_dynamic_frictions == old_dynamic_friction);
    CHECK(f.state.pairs.combined_static_frictions == old_static_friction);
    CHECK(f.state.pairs.combined_damping_ratios == old_damping);
    CHECK(f.state.pairs.combined_contact_frequencies == old_frequency);
    CHECK(f.state.pairs.combined_max_penetration_correction_speeds == old_correction_speed);
}
