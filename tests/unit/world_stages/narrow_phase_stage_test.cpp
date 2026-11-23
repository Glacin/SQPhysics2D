#include "world/world_stages.h"

#include <catch_amalgamated.hpp>
#include <limits>

namespace {

sqp::u32 add_body(sqp::WorldState& state, sqp::IslandId island, sqp::Vec2 position,
                  sqp::Vec2 velocity, sqp::BodyConfigs::CcdOverride mode, bool polygon = false,
                  sqp::real rotation = 0, sqp::real angular_velocity = 0, bool unit_box = false) {
    sqp::RigidBodyDef definition;
    definition.position = position;
    definition.rotation = rotation;
    definition.linear_velocity = velocity;
    definition.angular_velocity = angular_velocity;
    definition.configs.ccd_override = mode;
    const sqp::BodyHandle body = state.body_pool.add_rigid(definition);
    sqp::ShapeHandle shape;
    if (polygon) {
        sqp::PolygonShapeDef box;
        box.vertices = unit_box ? sqp::PolygonVertices{{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}
                                : sqp::PolygonVertices{{-2, -sqp::real(0.1)},
                                                       {2, -sqp::real(0.1)},
                                                       {2, sqp::real(0.1)},
                                                       {-2, sqp::real(0.1)}};
        shape = state.shape_pool.add_polygon(box);
    } else {
        sqp::CircleShapeDef circle;
        circle.radius = 1;
        shape = state.shape_pool.add_circle(circle);
    }
    REQUIRE(state.body_pool.attach_shape(body, shape, state.shape_pool));
    REQUIRE(state.island_pool.add_body(island, body, state.body_pool));
    const auto index = static_cast<sqp::u32>(state.broad_phase.bodies.size());
    state.broad_phase.bodies.push_back(body);
    state.broad_phase.shapes.push_back(shape);
    return index;
}

sqp::u32 add_lite_body(sqp::WorldState& state, sqp::IslandId island, sqp::Vec2 position,
                       sqp::Vec2 velocity, bool box) {
    sqp::ShapeHandle shape;
    sqp::BodyHandle body;
    if (box) {
        sqp::PolygonShapeDef geometry;
        geometry.vertices = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
        shape = state.shape_pool.add_polygon(geometry);
        sqp::AabbBoxBodyDef definition;
        definition.position = position;
        definition.linear_velocity = velocity;
        body = state.body_pool.add_aabb_box(definition, shape, state.shape_pool);
    } else {
        sqp::CircleShapeDef geometry;
        geometry.radius = 1;
        shape = state.shape_pool.add_circle(geometry);
        sqp::LiteCircleBodyDef definition;
        definition.position = position;
        definition.linear_velocity = velocity;
        body = state.body_pool.add_lite_circle(definition, shape, state.shape_pool);
    }
    REQUIRE(sqp::is_valid(body));
    REQUIRE(state.island_pool.add_body(island, body, state.body_pool));
    const auto index = static_cast<sqp::u32>(state.broad_phase.bodies.size());
    state.broad_phase.bodies.push_back(body);
    state.broad_phase.shapes.push_back(shape);
    return index;
}

void pair(sqp::WorldState& state, sqp::u32 a, sqp::u32 b) {
    state.broad_phase.candidates.push_back({a, b});
}

void clear_outputs(sqp::WorldState& state) {
    auto& output = state.narrow_phase;
    state.pairs.pairs->clear();
    state.pairs.manifolds->clear();
    output.toi_pairs.clear();
    output.toi_times.clear();
    output.toi_normals.clear();
    output.toi_impacts.clear();
}

void check_manifold_layout(const sqp::WorldState& state) {
    const auto& manifolds = *state.pairs.manifolds;
    CHECK(manifolds.size() == state.pairs.pairs->size());
    CHECK(manifolds.keys.size() == manifolds.counts.size() * 2);
    CHECK(manifolds.points_x.size() == manifolds.keys.size());
    CHECK(manifolds.points_y.size() == manifolds.keys.size());
    CHECK(manifolds.depths.size() == manifolds.keys.size());
    CHECK(manifolds.normals_x.size() == manifolds.counts.size());
    CHECK(manifolds.normals_y.size() == manifolds.counts.size());
}

} // namespace

TEST_CASE("TOI time rounds upward within the interval", "[world][narrow_phase][toi]") {
    CHECK(sqp::ceil_toi_time(10, 18, sqp::real(0.25)) == 12);
    CHECK(sqp::ceil_toi_time(10, 18, sqp::real(0.3)) == 13);
    CHECK(sqp::ceil_toi_time(10, 18, sqp::real(0)) == 10);
    CHECK(sqp::ceil_toi_time(10, 18, sqp::real(1)) == 18);
    for (const sqp::real fraction : {sqp::real(0), sqp::real(0.25), sqp::real(0.3), sqp::real(1)}) {
        const auto time = sqp::ceil_toi_time(10, 18, fraction);
        CHECK(time >= 10);
        CHECK(time <= 18);
        CHECK(static_cast<long double>(time - 10) >= static_cast<long double>(fraction) * 8);
    }
}

TEST_CASE("AABB body CCD uses fast swept bounds for boxes and lite circles",
          "[world][narrow_phase][ccd]") {
    const bool circle = GENERATE(false, true);
    const bool reversed = GENERATE(false, true);
    sqp::WorldState state;
    const auto island = state.island_pool.create();
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto moving = add_lite_body(state, island, {0, 0}, {10, 0}, !circle);
    const auto target = add_lite_body(state, island, {6, 0}, {}, true);
    pair(state, reversed ? target : moving, reversed ? moving : target);
    REQUIRE(sqp::run_narrow_phase(state));
    REQUIRE(state.narrow_phase.toi_pairs.size() == 1);
    CHECK(state.narrow_phase.toi_times[0] ==
          Catch::Approx(sqp::real(0.4) * sqp::real(sqp::config::TU_PER_SECOND)).margin(2));
    CHECK(state.narrow_phase.toi_normals[0].x == Catch::Approx(reversed ? 1 : -1));
}

TEST_CASE("contact manifolds use two point slots per pair", "[world][narrow_phase][manifold]") {
    sqp::ContactManifolds manifolds;
    sqp::ContactManifold single;
    single.count = 1;
    single.normal = {1, 0};
    single.points[0].point = {2, 3};
    single.points[0].depth = sqp::real(0.5);
    single.points[0].key.feature_a.type = sqp::ContactFeatureType::Circle;
    REQUIRE(manifolds.append(single));
    REQUIRE(manifolds.size() == 1);
    CHECK(manifolds.counts[0] == 1);
    CHECK(manifolds.points_x[0] == 2);
    CHECK(manifolds.points_y[0] == 3);
    CHECK(manifolds.depths[0] == sqp::real(0.5));
    CHECK(manifolds.keys[0] == single.points[0].key);
    CHECK(manifolds.points_x[1] == 0);
    CHECK(manifolds.points_y[1] == 0);
    CHECK(manifolds.depths[1] == 0);
    CHECK(manifolds.keys[1] == sqp::ContactPointKey{});
    CHECK(manifolds.at(0).points[0].key == single.points[0].key);

    sqp::ContactManifold double_point;
    double_point.count = 2;
    double_point.normal = {0, -1};
    double_point.points[0].point = {4, 5};
    double_point.points[1].point = {6, 7};
    double_point.points[1].depth = 2;
    double_point.points[1].key.feature_b.type = sqp::ContactFeatureType::Vertex;
    REQUIRE(manifolds.append(double_point));
    REQUIRE(manifolds.size() == 2);
    CHECK(manifolds.counts[1] == 2);
    CHECK(manifolds.points_x[2] == 4);
    CHECK(manifolds.points_x[3] == 6);
    CHECK(manifolds.points_y[3] == 7);
    CHECK(manifolds.depths[3] == 2);
    CHECK(manifolds.keys[3] == double_point.points[1].key);
    CHECK(manifolds.normals_x[1] == 0);
    CHECK(manifolds.normals_y[1] == -1);
    CHECK(manifolds.at(1).points[1].point == double_point.points[1].point);
    single.count = 0;
    CHECK_FALSE(manifolds.append(single));
    single.count = 3;
    CHECK_FALSE(manifolds.append(single));
    CHECK(manifolds.size() == 2);
    const auto key_capacity = manifolds.keys.capacity();
    manifolds.clear();
    CHECK(manifolds.empty());
    CHECK(manifolds.keys.empty());
    CHECK(manifolds.normals_x.empty());
    CHECK(manifolds.keys.capacity() == key_capacity);
}

TEST_CASE("narrow phase appends discrete pairs in candidate order", "[world][narrow_phase]") {
    sqp::WorldState state;
    const auto island = state.island_pool.create();
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto a = add_body(state, island, {0, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    const auto b = add_body(state, island, {1, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    const auto c = add_body(state, island, {8, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    auto versions = state.body_pool.motion_version_column(sqp::BodyConfigs::Kind::Rigid);
    versions[*state.body_pool.dense_index(state.broad_phase.bodies[a])] = sqp::u64{1} << 40;
    versions[*state.body_pool.dense_index(state.broad_phase.bodies[b])] = 73;
    pair(state, a, c);
    pair(state, a, b);
    pair(state, b, c);
    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.broad_phase.candidates.empty());
    check_manifold_layout(state);
    REQUIRE(state.pairs.pairs->size() == 1);
    REQUIRE(state.pairs.manifolds->size() == state.pairs.pairs->size());
    CHECK(state.pairs.manifolds->at(0).count == 1);
    CHECK(state.pairs.manifolds->at(0).normal.x == Catch::Approx(1));
    CHECK((*state.pairs.pairs)[0].body_a == state.broad_phase.bodies[a]);
    CHECK((*state.pairs.pairs)[0].body_b == state.broad_phase.bodies[b]);
    CHECK((*state.pairs.pairs)[0].motion_version_a == (sqp::u64{1} << 40));
    CHECK((*state.pairs.pairs)[0].motion_version_b == 73);
    CHECK(state.narrow_phase.toi_pairs.empty());
    versions[*state.body_pool.dense_index(state.broad_phase.bodies[b])] = 99;
    pair(state, a, b);
    REQUIRE(sqp::run_narrow_phase(state));
    check_manifold_layout(state);
    REQUIRE(state.pairs.pairs->size() == 2);
    REQUIRE(state.pairs.manifolds->size() == 2);
    CHECK((*state.pairs.pairs)[0].motion_version_b == 73);
    CHECK((*state.pairs.pairs)[1].motion_version_b == 99);
    CHECK((*state.pairs.pairs)[1].body_a == state.broad_phase.bodies[a]);
    CHECK((*state.pairs.pairs)[1].body_b == state.broad_phase.bodies[b]);
    CHECK(state.pairs.manifolds->at(1).count == state.pairs.manifolds->at(0).count);
    CHECK(state.pairs.manifolds->at(1).normal == state.pairs.manifolds->at(0).normal);
}

TEST_CASE("narrow phase appends contact manifolds and TOI records across runs",
          "[world][narrow_phase]") {
    sqp::WorldState state;
    const auto island = state.island_pool.create();
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto contact_a =
        add_body(state, island, {0, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    const auto contact_b =
        add_body(state, island, {1, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    const auto moving =
        add_body(state, island, {-10, 0}, {10, 0}, sqp::BodyConfigs::CcdOverride::Shapecast);
    const auto target =
        add_body(state, island, {0, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    pair(state, contact_a, contact_b);
    pair(state, moving, target);

    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.broad_phase.candidates.empty());
    check_manifold_layout(state);
    REQUIRE(state.pairs.pairs->size() == 1);
    REQUIRE(state.narrow_phase.toi_pairs.size() == 1);
    const auto first_contact = state.pairs.manifolds->at(0);
    const auto first_toi_time = state.narrow_phase.toi_times[0];
    const auto first_toi_normal = state.narrow_phase.toi_normals[0];
    auto versions = state.body_pool.motion_version_column(sqp::BodyConfigs::Kind::Rigid);
    versions[*state.body_pool.dense_index(state.broad_phase.bodies[target])] = 42;

    pair(state, contact_a, contact_b);
    pair(state, moving, target);
    REQUIRE(sqp::run_narrow_phase(state));
    check_manifold_layout(state);
    REQUIRE(state.pairs.pairs->size() == 2);
    REQUIRE(state.narrow_phase.toi_pairs.size() == 2);
    REQUIRE(state.narrow_phase.toi_times.size() == 2);
    REQUIRE(state.narrow_phase.toi_normals.size() == 2);

    CHECK((*state.pairs.pairs)[0].body_a == (*state.pairs.pairs)[1].body_a);
    CHECK((*state.pairs.pairs)[0].body_b == (*state.pairs.pairs)[1].body_b);
    CHECK(state.pairs.manifolds->at(0).count == first_contact.count);
    CHECK(state.pairs.manifolds->at(0).normal == first_contact.normal);
    CHECK(state.pairs.manifolds->at(0).points[0].key == first_contact.points[0].key);
    CHECK(state.pairs.manifolds->at(1).normal == first_contact.normal);
    CHECK(state.narrow_phase.toi_pairs[0].motion_version_b != 42);
    CHECK(state.narrow_phase.toi_pairs[1].motion_version_b == 42);
    CHECK(state.narrow_phase.toi_times[0] == first_toi_time);
    CHECK(state.narrow_phase.toi_times[1] == first_toi_time);
    CHECK(state.narrow_phase.toi_normals[0] == first_toi_normal);
    CHECK(state.narrow_phase.toi_normals[1] == first_toi_normal);

    CHECK(state.broad_phase.candidates.empty());
    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.pairs.pairs->size() == 2);
    CHECK(state.pairs.manifolds->size() == 2);
    CHECK(state.narrow_phase.toi_pairs.size() == 2);
    CHECK(state.narrow_phase.toi_times.size() == 2);
    CHECK(state.narrow_phase.toi_normals.size() == 2);
}

TEST_CASE("narrow phase retains candidates on failure and consumes them after retry",
          "[world][narrow_phase]") {
    sqp::WorldState state;
    const auto island = state.island_pool.create();
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto a = add_body(state, island, {0, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    const auto b = add_body(state, island, {1, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    pair(state, a, b);
    REQUIRE(sqp::run_narrow_phase(state));
    REQUIRE(state.pairs.pairs->size() == 1);
    CHECK(state.broad_phase.candidates.empty());

    pair(state, a, b);
    const auto shape = state.broad_phase.shapes.back();
    state.broad_phase.shapes.pop_back();
    CHECK_FALSE(sqp::run_narrow_phase(state));
    CHECK(state.broad_phase.candidates == std::vector<sqp::BroadPhase::Pair>{{a, b}});
    CHECK(state.pairs.pairs->size() == 1);
    CHECK(state.pairs.manifolds->size() == 1);

    state.broad_phase.shapes.push_back(shape);
    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.broad_phase.candidates.empty());
    CHECK(state.pairs.pairs->size() == 2);
    check_manifold_layout(state);
}

TEST_CASE("narrow phase records absolute TOI with B-to-A normal", "[world][narrow_phase]") {
    sqp::WorldState state;
    const auto early = state.island_pool.create();
    const auto late = state.island_pool.create();
    REQUIRE(state.island_pool.set_time(late, sqp::config::TU_PER_SECOND));
    state.interval.target_end = 2 * sqp::config::TU_PER_SECOND;
    const auto a =
        add_body(state, early, {-10, 0}, {5, 0}, sqp::BodyConfigs::CcdOverride::Shapecast);
    const auto b = add_body(state, late, {0, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    auto versions = state.body_pool.motion_version_column(sqp::BodyConfigs::Kind::Rigid);
    versions[*state.body_pool.dense_index(state.broad_phase.bodies[a])] = 55;
    versions[*state.body_pool.dense_index(state.broad_phase.bodies[b])] = sqp::u64{1} << 48;
    pair(state, a, b);
    REQUIRE(sqp::run_narrow_phase(state));
    check_manifold_layout(state);
    REQUIRE(state.narrow_phase.toi_pairs.size() == 1);
    REQUIRE(state.narrow_phase.toi_times.size() == 1);
    REQUIRE(state.narrow_phase.toi_normals.size() == 1);
    CHECK(state.narrow_phase.toi_pairs[0].shape_a == state.broad_phase.shapes[a]);
    CHECK(state.narrow_phase.toi_pairs[0].motion_version_a == 55);
    CHECK(state.narrow_phase.toi_pairs[0].motion_version_b == (sqp::u64{1} << 48));
    CHECK(state.narrow_phase.toi_times[0] ==
          Catch::Approx(sqp::real(1.6) * sqp::real(sqp::config::TU_PER_SECOND)).margin(2));
    CHECK(state.narrow_phase.toi_normals[0].x == Catch::Approx(-1));
    CHECK(state.pairs.pairs->empty());
    CHECK(state.pairs.manifolds->empty());
}

TEST_CASE("narrow phase reuses proxy data across pairs and refreshes it each run",
          "[world][narrow_phase]") {
    sqp::WorldState state;
    const auto early = state.island_pool.create();
    const auto late = state.island_pool.create();
    REQUIRE(state.island_pool.set_time(late, sqp::config::TU_PER_SECOND / 2));
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto moving =
        add_body(state, early, {-6, 0}, {8, 0}, sqp::BodyConfigs::CcdOverride::Shapecast);
    const auto first = add_body(state, early, {0, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    const auto second = add_body(state, late, {2, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    const auto invalid =
        add_body(state, early, {20, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    state.broad_phase.shapes[invalid] = state.broad_phase.shapes[first];
    const auto enqueue_candidates = [&] {
        pair(state, moving, second);
        pair(state, moving, invalid);
        pair(state, second, invalid);
        pair(state, moving, first);
    };
    enqueue_candidates();

    REQUIRE(sqp::run_narrow_phase(state));
    REQUIRE(state.narrow_phase.toi_pairs.size() == 2);
    REQUIRE(state.narrow_phase.toi_times.size() == 2);
    CHECK(state.narrow_phase.toi_pairs[0].body_b == state.broad_phase.bodies[second]);
    CHECK(state.narrow_phase.toi_pairs[1].body_b == state.broad_phase.bodies[first]);
    CHECK(state.narrow_phase.toi_times[0] ==
          Catch::Approx(sqp::real(0.75) * sqp::real(sqp::config::TU_PER_SECOND)).margin(2));
    CHECK(state.narrow_phase.toi_times[1] ==
          Catch::Approx(sqp::real(0.5) * sqp::real(sqp::config::TU_PER_SECOND)).margin(2));
    CHECK(state.narrow_phase.candidate_cache[invalid].epoch ==
          state.narrow_phase.candidate_cache_epoch);
    CHECK_FALSE(state.narrow_phase.candidate_cache[invalid].valid);
    const auto capacity = state.narrow_phase.candidate_cache.capacity();
    const auto first_epoch = state.narrow_phase.candidate_cache_epoch;

    REQUIRE(state.body_pool.set_velocity(state.broad_phase.bodies[moving], {}, 0));
    clear_outputs(state);
    enqueue_candidates();
    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.narrow_phase.candidate_cache.capacity() == capacity);
    CHECK(state.narrow_phase.candidate_cache_epoch != first_epoch);
    CHECK(state.narrow_phase.toi_pairs.empty());
    CHECK(state.pairs.pairs->empty());

    REQUIRE(state.body_pool.set_velocity(state.broad_phase.bodies[moving], {8, 0}, 0));
    REQUIRE(state.body_pool.set_ccd_override(state.broad_phase.bodies[moving],
                                             sqp::BodyConfigs::CcdOverride::Disabled));
    clear_outputs(state);
    enqueue_candidates();
    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.narrow_phase.toi_pairs.empty());
    REQUIRE(state.pairs.pairs->size() == 2);
    CHECK((*state.pairs.pairs)[0].body_b == state.broad_phase.bodies[second]);
    CHECK((*state.pairs.pairs)[1].body_b == state.broad_phase.bodies[first]);
    check_manifold_layout(state);

    state.narrow_phase.candidate_cache_epoch = std::numeric_limits<sqp::u32>::max();
    enqueue_candidates();
    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.narrow_phase.candidate_cache_epoch == 1);
    CHECK(state.narrow_phase.candidate_cache[moving].epoch == 1);
    CHECK(state.narrow_phase.candidate_cache[invalid].epoch == 1);
}

TEST_CASE("narrow phase reads motion versions from each body kind", "[world][narrow_phase]") {
    sqp::WorldState state;
    const auto island = state.island_pool.create();
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto rigid = add_body(state, island, {0, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    sqp::CircleShapeDef circle;
    circle.radius = 1;
    const auto shape = state.shape_pool.add_circle(circle);
    sqp::LiteCircleBodyDef definition;
    definition.position = {1, 0};
    const auto lite_body = state.body_pool.add_lite_circle(definition, shape, state.shape_pool);
    REQUIRE(state.island_pool.add_body(island, lite_body, state.body_pool));
    const auto lite = static_cast<sqp::u32>(state.broad_phase.bodies.size());
    state.broad_phase.bodies.push_back(lite_body);
    state.broad_phase.shapes.push_back(shape);
    state.body_pool.motion_version_column(sqp::BodyConfigs::Kind::Rigid)[0] = 7;
    state.body_pool.motion_version_column(sqp::BodyConfigs::Kind::LiteCircle)[0] = sqp::u64{1}
                                                                                   << 50;
    pair(state, rigid, lite);
    REQUIRE(sqp::run_narrow_phase(state));
    REQUIRE(state.pairs.pairs->size() == 1);
    CHECK(state.pairs.manifolds->size() == 1);
    CHECK((*state.pairs.pairs)[0].motion_version_a == 7);
    CHECK((*state.pairs.pairs)[0].motion_version_b == (sqp::u64{1} << 50));
}

TEST_CASE("narrow phase pairs identify attached shapes", "[world][narrow_phase]") {
    sqp::WorldState state;
    const auto island = state.island_pool.create();
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto a = add_body(state, island, {0, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    const auto b = add_body(state, island, {4, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    sqp::CircleShapeDef offset_circle;
    offset_circle.radius = 1;
    offset_circle.local_offset = {3, 0};
    const auto second_shape = state.shape_pool.add_circle(offset_circle);
    REQUIRE(
        state.body_pool.attach_shape(state.broad_phase.bodies[a], second_shape, state.shape_pool));
    state.broad_phase.shapes[a] = second_shape;
    pair(state, a, b);
    REQUIRE(sqp::run_narrow_phase(state));
    check_manifold_layout(state);
    REQUIRE(state.pairs.pairs->size() == 1);
    REQUIRE(state.pairs.manifolds->size() == 1);
    REQUIRE(state.pairs.manifolds->at(0).count == 1);
    CHECK((*state.pairs.pairs)[0].shape_a == second_shape);
    CHECK((*state.pairs.pairs)[0].shape_b == state.broad_phase.shapes[b]);
    CHECK(state.pairs.manifolds->at(0).points[0].key.feature_a.type ==
          sqp::ContactFeatureType::Circle);
    state.broad_phase.shapes[a] = state.broad_phase.shapes[b];
    clear_outputs(state);
    pair(state, a, b);
    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.pairs.pairs->empty());
    CHECK(state.pairs.manifolds->empty());
}

TEST_CASE("narrow phase builds polygon contact features", "[world][narrow_phase][manifold]") {
    sqp::WorldState state;
    const auto island = state.island_pool.create();
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto a = add_body(state, island, {0, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled,
                            true, 0, 0, true);
    const auto b = add_body(state, island, {sqp::real(1.5), 0}, {},
                            sqp::BodyConfigs::CcdOverride::Disabled, true, 0, 0, true);
    pair(state, a, b);
    REQUIRE(sqp::run_narrow_phase(state));
    REQUIRE(state.pairs.pairs->size() == 1);
    REQUIRE(state.pairs.manifolds->size() == 1);
    const auto manifold = state.pairs.manifolds->at(0);
    REQUIRE(manifold.count == 2);
    CHECK(manifold.normal.x == Catch::Approx(1));
    for (const auto& contact : manifold.points) {
        CHECK(contact.key.feature_a.type != sqp::ContactFeatureType::Circle);
        CHECK(contact.key.feature_b.type != sqp::ContactFeatureType::Circle);
    }
}

TEST_CASE("narrow phase builds circle polygon contacts in both orders",
          "[world][narrow_phase][manifold]") {
    sqp::WorldState state;
    const auto island = state.island_pool.create();
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto polygon = add_body(state, island, {0, 0}, {},
                                  sqp::BodyConfigs::CcdOverride::Disabled, true, 0, 0, true);
    const auto circle =
        add_body(state, island, {sqp::real(1.5), 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    pair(state, polygon, circle);
    pair(state, circle, polygon);
    REQUIRE(sqp::run_narrow_phase(state));
    check_manifold_layout(state);
    REQUIRE(state.pairs.pairs->size() == 2);
    REQUIRE(state.pairs.manifolds->size() == 2);
    const auto forward = state.pairs.manifolds->at(0);
    const auto reverse = state.pairs.manifolds->at(1);
    REQUIRE(forward.count == 1);
    REQUIRE(reverse.count == 1);
    CHECK(forward.normal.x == Catch::Approx(1));
    CHECK(reverse.normal.x == Catch::Approx(-1));
    CHECK(forward.points[0].key.feature_a.type == sqp::ContactFeatureType::Edge);
    CHECK(forward.points[0].key.feature_b.type == sqp::ContactFeatureType::Circle);
    CHECK(reverse.points[0].key.feature_a.type == sqp::ContactFeatureType::Circle);
    CHECK(reverse.points[0].key.feature_b.type == sqp::ContactFeatureType::Edge);
}

TEST_CASE("Auto thresholds and family switches route to discrete or TOI", "[world][narrow_phase]") {
    sqp::WorldState state;
    const auto island = state.island_pool.create();
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto a = add_body(state, island, {-sqp::real(2.2), 0}, {sqp::real(0.39), 0},
                            sqp::BodyConfigs::CcdOverride::Auto);
    const auto b = add_body(state, island, {0, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    pair(state, a, b);
    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.pairs.pairs->size() == 1);
    CHECK(state.pairs.manifolds->size() == 1);
    CHECK(state.narrow_phase.toi_pairs.empty());

    REQUIRE(state.body_pool.set_velocity(state.broad_phase.bodies[a], {sqp::real(0.41), 0}, 0));
    clear_outputs(state);
    pair(state, a, b);
    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.pairs.pairs->empty());
    CHECK(state.pairs.manifolds->empty());
    CHECK(state.narrow_phase.toi_pairs.size() == 1);

    state.settings.ccd.enable_circle = false;
    clear_outputs(state);
    pair(state, a, b);
    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.pairs.pairs->size() == 1);
    CHECK(state.pairs.manifolds->size() == 1);
    CHECK(state.narrow_phase.toi_pairs.empty());

    state.settings.ccd.enable_circle = true;
    REQUIRE(state.body_pool.set_velocity(state.broad_phase.bodies[a], {sqp::real(0.39), 0}, 0));
    clear_outputs(state);
    pair(state, a, b);
    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.pairs.pairs->size() == 1);
    CHECK(state.narrow_phase.toi_pairs.empty());
    REQUIRE(state.body_pool.set_velocity(state.broad_phase.bodies[a],
                                         {sqp::real(0.3), sqp::real(0.3)}, 0));
    clear_outputs(state);
    pair(state, a, b);
    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.pairs.pairs->empty());
    CHECK(state.narrow_phase.toi_pairs.size() == 1);
}

TEST_CASE("separating initial overlap records a later rotational reentry",
          "[world][narrow_phase]") {
    sqp::WorldState state;
    const auto island = state.island_pool.create();
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto blade = add_body(state, island, {0, 0}, {}, sqp::BodyConfigs::CcdOverride::Shapecast,
                                true, sqp::math::PI * sqp::real(0.5), sqp::math::PI * sqp::real(2));
    // Offset the target so the initial separating direction is not symmetric.
    const auto target = add_body(state, island, {sqp::real(0.5), sqp::real(1.5)}, {},
                                 sqp::BodyConfigs::CcdOverride::Disabled, true, 0, 0, true);
    pair(state, blade, target);
    REQUIRE(sqp::run_narrow_phase(state));
    REQUIRE(state.narrow_phase.toi_times.size() == 1);
    CHECK(state.narrow_phase.toi_times[0] > sqp::config::TU_PER_SECOND / 4);
    CHECK(state.narrow_phase.toi_times[0] < 3 * sqp::config::TU_PER_SECOND / 4);
}

TEST_CASE("rotational Auto and forced CCD use ShapeCast", "[world][narrow_phase]") {
    sqp::WorldState state;
    const auto island = state.island_pool.create();
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto blade = add_body(state, island, {0, 0}, {}, sqp::BodyConfigs::CcdOverride::Auto,
                                true, 0, sqp::math::PI * sqp::real(0.5));
    const auto target = add_body(state, island, {0, sqp::real(1.5)}, {},
                                 sqp::BodyConfigs::CcdOverride::Disabled, true);
    pair(state, blade, target);
    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.narrow_phase.toi_pairs.size() == 1);
    CHECK(state.pairs.manifolds->empty());

    // The blade's roughly 2 m sweep radius turns a half-turn into over 2 m of
    // bounded motion, so the shared distance threshold still selects shapecast.
    state.settings.ccd.linear_start_threshold_m = sqp::real(2);
    clear_outputs(state);
    pair(state, blade, target);
    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.narrow_phase.toi_pairs.size() == 1);

    state.settings.ccd.linear_start_threshold_m = sqp::real(4);
    clear_outputs(state);
    pair(state, blade, target);
    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.narrow_phase.toi_pairs.empty());
    CHECK(state.pairs.pairs->size() == 1);
    REQUIRE(state.pairs.manifolds->size() == 1);
    CHECK(state.pairs.manifolds->at(0).count > 0);
    state.settings.ccd.linear_start_threshold_m = sqp::real(0.4);

    state.settings.ccd.enable_polygon = false;
    clear_outputs(state);
    pair(state, blade, target);
    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.narrow_phase.toi_pairs.empty());
    state.settings.ccd.enable_polygon = true;

    clear_outputs(state);
    pair(state, blade, target);
    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.narrow_phase.toi_pairs.size() == 1);

    REQUIRE(state.body_pool.set_ccd_override(state.broad_phase.bodies[target],
                                             sqp::BodyConfigs::CcdOverride::Shapecast));
    clear_outputs(state);
    pair(state, blade, target);
    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.narrow_phase.toi_pairs.size() == 1);
}

TEST_CASE("initial overlap keeps approaching TOI and rejects separating TOI",
          "[world][narrow_phase]") {
    sqp::WorldState state;
    const auto island = state.island_pool.create();
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto a =
        add_body(state, island, {-1, 0}, {1, 0}, sqp::BodyConfigs::CcdOverride::Shapecast);
    const auto b = add_body(state, island, {0, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    pair(state, a, b);
    REQUIRE(sqp::run_narrow_phase(state));
    REQUIRE(state.narrow_phase.toi_times.size() == 1);
    CHECK(state.narrow_phase.toi_times[0] == 0);
    REQUIRE(state.body_pool.set_velocity(state.broad_phase.bodies[a], {-4, 0}, 0));
    clear_outputs(state);
    pair(state, a, b);
    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.narrow_phase.toi_pairs.empty());
}

TEST_CASE("narrow phase ignores malformed candidates and still consumes the batch",
          "[world][narrow_phase]") {
    sqp::WorldState state;
    const auto island = state.island_pool.create();
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto a = add_body(state, island, {0, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    const auto b = add_body(state, island, {1, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    const auto alias = static_cast<sqp::u32>(state.broad_phase.bodies.size());
    state.broad_phase.bodies.push_back(state.broad_phase.bodies[a]);
    state.broad_phase.shapes.push_back(state.broad_phase.shapes[a]);
    const auto out_of_range = static_cast<sqp::u32>(state.broad_phase.bodies.size());

    pair(state, a, a);
    pair(state, a, alias);
    pair(state, out_of_range, b);
    pair(state, a, out_of_range);
    pair(state, a, b);
    REQUIRE(sqp::run_narrow_phase(state));
    CHECK(state.broad_phase.candidates.empty());
    REQUIRE(state.pairs.pairs->size() == 1);
    CHECK((*state.pairs.pairs)[0].body_a == state.broad_phase.bodies[a]);
    CHECK((*state.pairs.pairs)[0].body_b == state.broad_phase.bodies[b]);
    CHECK(state.narrow_phase.toi_pairs.empty());
    check_manifold_layout(state);
}

TEST_CASE("discrete narrow phase predicts each body from its own island time",
          "[world][narrow_phase]") {
    sqp::WorldState state;
    const auto early = state.island_pool.create();
    const auto late = state.island_pool.create();
    const auto half_second = sqp::config::TU_PER_SECOND / 2;
    REQUIRE(state.island_pool.set_time(late, half_second));
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto a = add_body(state, early, {0, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    const auto b = add_body(state, late, {-4, 0}, {10, 0}, sqp::BodyConfigs::CcdOverride::Disabled);
    pair(state, a, b);
    REQUIRE(sqp::run_narrow_phase(state));
    REQUIRE(state.pairs.pairs->size() == 1);
    CHECK(state.pairs.manifolds->size() == 1);
    CHECK(state.pairs.manifolds->at(0).normal.x == Catch::Approx(1));
    CHECK(state.narrow_phase.toi_pairs.empty());
}

TEST_CASE("reversing a CCD candidate reverses its TOI normal", "[world][narrow_phase]") {
    sqp::WorldState state;
    const auto island = state.island_pool.create();
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto moving =
        add_body(state, island, {-4, 0}, {4, 0}, sqp::BodyConfigs::CcdOverride::Shapecast);
    const auto target =
        add_body(state, island, {0, 0}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    pair(state, moving, target);
    pair(state, target, moving);
    REQUIRE(sqp::run_narrow_phase(state));
    REQUIRE(state.narrow_phase.toi_pairs.size() == 2);
    REQUIRE(state.narrow_phase.toi_times.size() == 2);
    REQUIRE(state.narrow_phase.toi_normals.size() == 2);
    CHECK(state.narrow_phase.toi_pairs[0].body_a == state.broad_phase.bodies[moving]);
    CHECK(state.narrow_phase.toi_pairs[1].body_a == state.broad_phase.bodies[target]);
    CHECK(state.narrow_phase.toi_times[0] == state.narrow_phase.toi_times[1]);
    CHECK(state.narrow_phase.toi_times[0] ==
          Catch::Approx(sqp::real(0.5) * sqp::real(sqp::config::TU_PER_SECOND)).margin(2));
    CHECK(state.narrow_phase.toi_normals[0].x == Catch::Approx(-1));
    CHECK(state.narrow_phase.toi_normals[1].x == Catch::Approx(1));
    CHECK(state.pairs.pairs->empty());
}

TEST_CASE("straight circle centers produce CCD for lite and rigid bodies",
          "[world][narrow_phase][ccd]") {
    for (const auto [rigid_a, rigid_b] : {std::pair{false, false}, std::pair{false, true},
                                          std::pair{true, false}, std::pair{true, true}}) {
        sqp::WorldState state;
        const auto island = state.island_pool.create();
        state.interval.target_end = sqp::config::TU_PER_SECOND;
        const auto add_circle = [&](bool rigid, sqp::Vec2 position, sqp::Vec2 velocity) {
            return rigid
                       ? add_body(state, island, position, velocity,
                                  sqp::BodyConfigs::CcdOverride::Shapecast, false, 0, sqp::real(1))
                       : add_lite_body(state, island, position, velocity, false);
        };
        const auto moving = add_circle(rigid_a, {-4, 0}, {4, 0});
        const auto target = add_circle(rigid_b, {0, 0}, {});
        pair(state, moving, target);
        pair(state, target, moving);
        REQUIRE(sqp::run_narrow_phase(state));
        REQUIRE(state.narrow_phase.toi_pairs.size() == 2);
        CHECK(state.narrow_phase.toi_times[0] == state.narrow_phase.toi_times[1]);
        CHECK(state.narrow_phase.toi_times[0] ==
              Catch::Approx(sqp::real(0.5) * sqp::real(sqp::config::TU_PER_SECOND)).margin(2));
        CHECK(state.narrow_phase.toi_normals[0].x == Catch::Approx(-1));
        CHECK(state.narrow_phase.toi_normals[1].x == Catch::Approx(1));
    }
}

TEST_CASE("LiteCircle sweeps against a translating polygon in both orders",
          "[world][narrow_phase][ccd]") {
    sqp::WorldState state;
    const auto island = state.island_pool.create();
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto circle = add_lite_body(state, island, {-4, 0}, {4, 0}, false);
    const auto polygon = add_body(state, island, {0, 0}, {1, 0},
                                  sqp::BodyConfigs::CcdOverride::Disabled, true, 0, 0, true);
    pair(state, circle, polygon);
    pair(state, polygon, circle);
    REQUIRE(sqp::run_narrow_phase(state));
    REQUIRE(state.narrow_phase.toi_pairs.size() == 2);
    CHECK(state.narrow_phase.toi_times[0] == state.narrow_phase.toi_times[1]);
    CHECK(state.narrow_phase.toi_times[0] ==
          Catch::Approx(sqp::real(2.0 / 3.0) * sqp::real(sqp::config::TU_PER_SECOND)).margin(2));
    CHECK(state.narrow_phase.toi_normals[0].x == Catch::Approx(-1));
    CHECK(state.narrow_phase.toi_normals[1].x == Catch::Approx(1));
}

TEST_CASE("rotating offset circles retain curved CCD motion", "[world][narrow_phase][ccd]") {
    sqp::WorldState state;
    const auto island = state.island_pool.create();
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto moving = add_body(state, island, {}, {}, sqp::BodyConfigs::CcdOverride::Shapecast,
                                 false, 0, sqp::math::PI * sqp::real(0.5));
    const sqp::real halfway = sqp::real(3) / sqp::math::sqrt(sqp::real(2));
    const auto target =
        add_body(state, island, {halfway, halfway}, {}, sqp::BodyConfigs::CcdOverride::Disabled);
    sqp::CircleShapeDef small_circle;
    small_circle.radius = sqp::real(0.2);
    small_circle.local_offset = {3, 0};
    const auto moving_shape = state.shape_pool.add_circle(small_circle);
    REQUIRE(state.body_pool.attach_shape(state.broad_phase.bodies[moving], moving_shape,
                                         state.shape_pool));
    state.broad_phase.shapes[moving] = moving_shape;
    small_circle.local_offset = {};
    const auto target_shape = state.shape_pool.add_circle(small_circle);
    REQUIRE(state.body_pool.attach_shape(state.broad_phase.bodies[target], target_shape,
                                         state.shape_pool));
    state.broad_phase.shapes[target] = target_shape;
    pair(state, moving, target);
    REQUIRE(sqp::run_narrow_phase(state));
    REQUIRE(state.narrow_phase.toi_pairs.size() == 1);
    CHECK(state.narrow_phase.toi_times[0] > sqp::config::TU_PER_SECOND / 3);
    CHECK(state.narrow_phase.toi_times[0] < 2 * sqp::config::TU_PER_SECOND / 3);
}

TEST_CASE("LiteCircle uses shape cast against a rotating polygon", "[world][narrow_phase][ccd]") {
    sqp::WorldState state;
    const auto island = state.island_pool.create();
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto blade = add_body(state, island, {}, {}, sqp::BodyConfigs::CcdOverride::Shapecast,
                                true, 0, sqp::math::PI * sqp::real(0.5));
    const auto circle = add_lite_body(state, island, {0, sqp::real(2.5)}, {}, false);
    pair(state, blade, circle);
    pair(state, circle, blade);
    REQUIRE(sqp::run_narrow_phase(state));
    REQUIRE(state.narrow_phase.toi_pairs.size() == 2);
    CHECK(state.narrow_phase.toi_times[0] == state.narrow_phase.toi_times[1]);
    CHECK(state.narrow_phase.toi_times[0] > 0);
    CHECK(state.narrow_phase.toi_times[0] < sqp::config::TU_PER_SECOND);
    CHECK(state.narrow_phase.toi_normals[0].y < 0);
    CHECK(state.narrow_phase.toi_normals[1].y > 0);
}

TEST_CASE("rotating rigid circle keeps its arc against a polygon", "[world][narrow_phase][ccd]") {
    sqp::WorldState state;
    const auto island = state.island_pool.create();
    state.interval.target_end = sqp::config::TU_PER_SECOND;
    const auto moving = add_body(state, island, {}, {}, sqp::BodyConfigs::CcdOverride::Shapecast,
                                 false, 0, sqp::math::PI * sqp::real(0.5));
    const sqp::real halfway = sqp::real(3) / sqp::math::sqrt(sqp::real(2));
    const auto target = add_body(state, island, {halfway, halfway}, {},
                                 sqp::BodyConfigs::CcdOverride::Disabled, true);
    sqp::CircleShapeDef circle;
    circle.radius = sqp::real(0.2);
    circle.local_offset = {3, 0};
    const auto moving_shape = state.shape_pool.add_circle(circle);
    REQUIRE(state.body_pool.attach_shape(state.broad_phase.bodies[moving], moving_shape,
                                         state.shape_pool));
    state.broad_phase.shapes[moving] = moving_shape;
    sqp::PolygonShapeDef polygon;
    polygon.vertices = {{-sqp::real(0.2), -sqp::real(0.2)},
                        {sqp::real(0.2), -sqp::real(0.2)},
                        {sqp::real(0.2), sqp::real(0.2)},
                        {-sqp::real(0.2), sqp::real(0.2)}};
    const auto target_shape = state.shape_pool.add_polygon(polygon);
    REQUIRE(state.body_pool.attach_shape(state.broad_phase.bodies[target], target_shape,
                                         state.shape_pool));
    state.broad_phase.shapes[target] = target_shape;
    pair(state, moving, target);
    REQUIRE(sqp::run_narrow_phase(state));
    REQUIRE(state.narrow_phase.toi_pairs.size() == 1);
    CHECK(state.narrow_phase.toi_times[0] > sqp::config::TU_PER_SECOND / 3);
    CHECK(state.narrow_phase.toi_times[0] < 2 * sqp::config::TU_PER_SECOND / 3);
}
