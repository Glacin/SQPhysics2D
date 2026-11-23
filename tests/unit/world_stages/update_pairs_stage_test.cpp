#include "world/world_stages.h"

#include <catch_amalgamated.hpp>

namespace {

struct Fixture {
    sqp::WorldState state;

    sqp::BodyHandle body() {
        return state.body_pool.add_rigid();
    }

    sqp::Pair append(sqp::BodyHandle a, sqp::BodyHandle b, sqp::real marker, sqp::i32 points = 1) {
        const auto versions_a = state.body_pool.motion_version_column(a.kind());
        const auto versions_b = state.body_pool.motion_version_column(b.kind());
        sqp::Pair pair{a,
                       b,
                       {},
                       {},
                       versions_a[*state.body_pool.dense_index(a)],
                       versions_b[*state.body_pool.dense_index(b)]};
        state.pairs.pairs->push_back(pair);
        sqp::ContactManifold manifold;
        manifold.count = points;
        manifold.normal = {marker, -marker};
        for (sqp::i32 slot = 0; slot < points; ++slot) {
            auto& point = manifold.points[slot];
            point.point = {marker + sqp::real(slot), marker * 2 + sqp::real(slot)};
            point.depth = marker + sqp::real(slot) / 10;
            point.key.feature_a.type = sqp::ContactFeatureType::Vertex;
            point.key.feature_a.index = static_cast<sqp::u16>(marker * 10 + slot);
            point.key.feature_b.type = sqp::ContactFeatureType::Edge;
            point.key.feature_b.index = static_cast<sqp::u16>(marker * 10 + slot + 1);
        }
        REQUIRE(state.pairs.manifolds->append(manifold));
        return pair;
    }

    void change_version(sqp::BodyHandle body) {
        auto versions = state.body_pool.motion_version_column(body.kind());
        ++versions[*state.body_pool.dense_index(body)];
    }
};

} // namespace

TEST_CASE("UpdatePairs accepts an empty ordinary pair set", "[world][update_pairs]") {
    Fixture f;
    REQUIRE(sqp::run_update_pairs(f.state));
    CHECK(f.state.pairs.pairs->empty());
    CHECK(f.state.pairs.manifolds->empty());
}

TEST_CASE("UpdatePairs removes a pair when either motion version changes",
          "[world][update_pairs]") {
    for (const bool change_a : {false, true}) {
        for (const bool change_b : {false, true}) {
            Fixture f;
            const auto a = f.body();
            const auto b = f.body();
            f.append(a, b, 1);
            if (change_a)
                f.change_version(a);
            if (change_b)
                f.change_version(b);
            REQUIRE(sqp::run_update_pairs(f.state));
            CHECK(f.state.pairs.pairs->size() == (change_a || change_b ? 0 : 1));
            CHECK(f.state.pairs.manifolds->size() == f.state.pairs.pairs->size());
        }
    }
}

TEST_CASE("UpdatePairs removes pairs whose body was deleted", "[world][update_pairs]") {
    for (const bool remove_a : {false, true}) {
        Fixture f;
        const auto a = f.body();
        const auto b = f.body();
        f.append(a, b, 2);
        REQUIRE(f.state.body_pool.remove(remove_a ? a : b));
        REQUIRE(sqp::run_update_pairs(f.state));
        CHECK(f.state.pairs.pairs->empty());
        CHECK(f.state.pairs.manifolds->empty());
    }
}

TEST_CASE("UpdatePairs preserves surviving order and every manifold column",
          "[world][update_pairs]") {
    Fixture f;
    const auto a = f.body();
    const auto b = f.body();
    const auto c = f.body();
    const auto d = f.body();
    const auto e = f.body();
    const auto first = f.append(a, b, 1, 2);
    f.append(a, c, 2);
    const auto middle = f.append(b, d, 3);
    f.append(c, e, 4, 2);
    const auto last = f.append(d, e, 5, 2);
    const auto first_manifold = f.state.pairs.manifolds->at(0);
    const auto middle_manifold = f.state.pairs.manifolds->at(2);
    const auto last_manifold = f.state.pairs.manifolds->at(4);
    const auto key_capacity = f.state.pairs.manifolds->keys.capacity();
    f.change_version(c);

    REQUIRE(sqp::run_update_pairs(f.state));
    const auto& pairs = *f.state.pairs.pairs;
    const auto& manifolds = *f.state.pairs.manifolds;
    REQUIRE(pairs.size() == 3);
    CHECK(pairs[0].body_a == first.body_a);
    CHECK(pairs[1].body_b == middle.body_b);
    CHECK(pairs[2].body_b == last.body_b);
    CHECK(manifolds.size() == pairs.size());
    CHECK(manifolds.keys.capacity() == key_capacity);
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        const auto actual = manifolds.at(i);
        const auto expected = i == 0 ? first_manifold : i == 1 ? middle_manifold : last_manifold;
        CHECK(actual.count == expected.count);
        CHECK(actual.normal == expected.normal);
        for (std::size_t slot = 0; slot < 2; ++slot) {
            CHECK(actual.points[slot].point == expected.points[slot].point);
            CHECK(actual.points[slot].depth == expected.points[slot].depth);
            CHECK(actual.points[slot].key == expected.points[slot].key);
        }
    }
    CHECK(manifolds.keys.size() == 2 * pairs.size());
}

TEST_CASE("UpdatePairs reads the version column for each body kind", "[world][update_pairs]") {
    Fixture f;
    const auto rigid = f.body();
    sqp::CircleShapeDef circle;
    circle.radius = 1;
    const auto shape = f.state.shape_pool.add_circle(circle);
    const auto lite = f.state.body_pool.add_lite_circle({}, shape, f.state.shape_pool);
    REQUIRE(f.state.body_pool.is_valid(lite));
    f.state.body_pool.motion_version_column(
        sqp::BodyConfigs::Kind::Rigid)[*f.state.body_pool.dense_index(rigid)] = 7;
    f.state.body_pool.motion_version_column(
        sqp::BodyConfigs::Kind::LiteCircle)[*f.state.body_pool.dense_index(lite)] = sqp::u64{1}
                                                                                    << 50;
    f.append(rigid, lite, 1);
    REQUIRE(sqp::run_update_pairs(f.state));
    CHECK(f.state.pairs.pairs->size() == 1);
    f.change_version(lite);
    REQUIRE(sqp::run_update_pairs(f.state));
    CHECK(f.state.pairs.pairs->empty());
    CHECK(f.state.pairs.manifolds->empty());
}

TEST_CASE("UpdatePairs rejects inconsistent manifold columns without filtering",
          "[world][update_pairs]") {
    for (int column = 0; column < 8; ++column) {
        Fixture f;
        f.append(f.body(), f.body(), 1);
        auto& m = *f.state.pairs.manifolds;
        switch (column) {
        case 0:
            m.counts.pop_back();
            break;
        case 1:
            m.normals_x.pop_back();
            break;
        case 2:
            m.normals_y.pop_back();
            break;
        case 3:
            m.keys.pop_back();
            break;
        case 4:
            m.points_x.pop_back();
            break;
        case 5:
            m.points_y.pop_back();
            break;
        case 6:
            m.depths.pop_back();
            break;
        default:
            f.state.pairs.pairs->pop_back();
            break;
        }
        const auto pair_count = f.state.pairs.pairs->size();
        CHECK_FALSE(sqp::run_update_pairs(f.state));
        CHECK(f.state.pairs.pairs->size() == pair_count);
    }
}

TEST_CASE("microstep interval bypasses UpdatePairs even with no contacts", "[world][update_pairs]") {
    Fixture f;
    CHECK(sqp::should_run_update_pairs(f.state));
    f.state.interval.microstep = true;
    CHECK_FALSE(sqp::should_run_update_pairs(f.state));
    f.state.interval.microstep = false;
    CHECK(sqp::should_run_update_pairs(f.state));
}
