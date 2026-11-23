#include "sqphysics2d/detail/data_storage.h"

#include <algorithm>
#include <catch_amalgamated.hpp>

namespace {

std::size_t count(const sqp::Island& island) {
    if (const auto* small = std::get_if<sqp::Island::SmallIsland>(&island.bodies))
        return static_cast<std::size_t>(
            std::find(small->data.begin(), small->data.end(), sqp::BodyHandle{}) -
            small->data.begin());
    return std::get<std::vector<sqp::BodyHandle>>(island.bodies).size();
}

} // namespace

TEST_CASE("Island TOI counts reset on tick and slot reuse", "[body][island]") {
    sqp::IslandPool islands;
    const auto first = islands.create();
    const auto second = islands.create();
    REQUIRE(islands.set_toi_truncations(first, 3));
    REQUIRE(islands.increment_toi_truncations(second));
    CHECK(islands.toi_truncations(first) == 3);
    CHECK(islands.toi_truncations(second) == 1);
    islands.reset_toi_truncations();
    CHECK(islands.toi_truncations(first) == 0);
    CHECK(islands.toi_truncations(second) == 0);
    REQUIRE(islands.set_toi_truncations(first, 4));
    REQUIRE(islands.destroy(first));
    const auto reused = islands.create();
    CHECK(reused.index == first.index);
    CHECK(reused.generation != first.generation);
    CHECK(islands.toi_truncations(reused) == 0);
    CHECK_FALSE(islands.increment_toi_truncations(first));
}

TEST_CASE("island membership follows all body kinds and changes representation", "[body][island]") {
    sqp::ShapePool shapes;
    const auto circle_shape = shapes.add_circle({1});
    sqp::PolygonShapeDef box_def;
    box_def.vertices = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
    const auto box_shape = shapes.add_polygon(box_def);
    sqp::BodyPool bodies;
    const auto rigid = bodies.add_rigid();
    const auto extra = bodies.add_rigid();
    const auto circle = bodies.add_lite_circle({}, circle_shape, shapes);
    const auto box = bodies.add_aabb_box({}, box_shape, shapes);
    sqp::IslandPool islands;
    const auto first = islands.create();
    const auto second = islands.create();
    REQUIRE(first.valid());
    REQUIRE(second.valid());
    CHECK_FALSE(bodies.island_of(rigid));
    CHECK_FALSE(bodies.island_of({}));
    REQUIRE(islands.add_body(first, rigid, bodies));
    REQUIRE(islands.add_body(first, circle, bodies));
    REQUIRE(islands.add_body(first, circle, bodies));
    REQUIRE(islands.get(first));
    CHECK(std::holds_alternative<sqp::Island::SmallIsland>(islands.get(first)->bodies));
    CHECK(count(*islands.get(first)) == 2);
    REQUIRE(islands.add_body(first, box, bodies));
    CHECK(std::holds_alternative<sqp::Island::SmallIsland>(islands.get(first)->bodies));
    REQUIRE(islands.add_body(first, extra, bodies));
    CHECK(std::holds_alternative<std::vector<sqp::BodyHandle>>(islands.get(first)->bodies));
    CHECK(count(*islands.get(first)) == 4);
    CHECK(bodies.island_of(rigid) == first);
    CHECK(bodies.island_of(circle) == first);
    CHECK(bodies.island_of(box) == first);

    REQUIRE(islands.add_body(second, box, bodies));
    CHECK(count(*islands.get(first)) == 3);
    CHECK(count(*islands.get(second)) == 1);
    CHECK(bodies.island_of(box) == second);
    REQUIRE(islands.remove_body(first, circle, bodies));
    CHECK(std::holds_alternative<sqp::Island::SmallIsland>(islands.get(first)->bodies));
    CHECK_FALSE(bodies.island_of(circle));
    CHECK_FALSE(islands.remove_body(first, circle, bodies));
    REQUIRE(islands.set_time(second, 42));
    CHECK(islands.get(second)->time == 42);
}

TEST_CASE("linked body and island removal reject stale handles", "[body][island]") {
    sqp::BodyPool bodies;
    sqp::IslandPool islands;
    const auto body = bodies.add_rigid();
    const auto island = islands.create();
    REQUIRE(islands.add_body(island, body, bodies));
    CHECK_FALSE(bodies.remove(body));
    CHECK_FALSE(islands.destroy(island));
    REQUIRE(bodies.remove(body, islands));
    CHECK_FALSE(bodies.island_of(body));
    CHECK_FALSE(islands.has_body(island, body));
    REQUIRE(islands.destroy(island));
    CHECK_FALSE(islands.get(island));
    const auto reused_body = bodies.add_rigid();
    const auto reused_island = islands.create();
    CHECK(reused_body.index() == body.index());
    CHECK(reused_body.generation != body.generation);
    CHECK(reused_island.index == island.index);
    CHECK(reused_island.generation != island.generation);
    CHECK_FALSE(islands.add_body(island, reused_body, bodies));
    CHECK_FALSE(islands.add_body(reused_island, body, bodies));
    CHECK_FALSE(bodies.island_of(reused_body));
}

TEST_CASE("destroying an island detaches all members", "[body][island]") {
    sqp::BodyPool bodies;
    sqp::IslandPool islands;
    const auto island = islands.create();
    const auto a = bodies.add_rigid();
    const auto b = bodies.add_rigid();
    const auto c = bodies.add_rigid();
    REQUIRE(islands.add_body(island, a, bodies));
    REQUIRE(islands.add_body(island, b, bodies));
    REQUIRE(islands.add_body(island, c, bodies));
    REQUIRE(islands.destroy(island, bodies));
    CHECK_FALSE(bodies.island_of(a));
    CHECK_FALSE(bodies.island_of(b));
    CHECK_FALSE(bodies.island_of(c));
    CHECK_FALSE(islands.get(island));
    CHECK(bodies.is_valid(a));
}

TEST_CASE("linked clear preserves islands and invalidates body handles", "[body][island]") {
    sqp::BodyPool bodies;
    sqp::IslandPool islands;
    const auto island = islands.create();
    const auto body = bodies.add_rigid();
    REQUIRE(islands.add_body(island, body, bodies));
    CHECK_FALSE(bodies.clear());
    CHECK(bodies.is_valid(body));
    CHECK(islands.has_body(island, body));
    REQUIRE(bodies.clear(islands));
    CHECK_FALSE(bodies.is_valid(body));
    CHECK(count(*islands.get(island)) == 0);
    const auto reused = bodies.add_rigid();
    CHECK(reused.index() == body.index());
    CHECK(reused.generation != body.generation);
    CHECK_FALSE(bodies.island_of(reused));
}

TEST_CASE("membership survives dense swaps and pool moves", "[body][island]") {
    sqp::BodyPool bodies;
    sqp::IslandPool islands;
    const auto island = islands.create();
    const auto first = bodies.add_rigid();
    const auto removed = bodies.add_rigid();
    const auto last = bodies.add_rigid();
    REQUIRE(islands.add_body(island, first, bodies));
    REQUIRE(islands.add_body(island, last, bodies));
    REQUIRE(bodies.remove(removed));
    CHECK(bodies.dense_index(last) == 1);
    sqp::BodyPool moved_bodies = std::move(bodies);
    sqp::IslandPool moved_islands = std::move(islands);
    CHECK(moved_bodies.island_of(last) == island);
    CHECK(moved_islands.has_body(island, last));
    REQUIRE(moved_bodies.remove(last, moved_islands));
    CHECK(moved_bodies.island_of(first) == island);
}

TEST_CASE("cross-pool mismatch leaves membership intact", "[body][island]") {
    sqp::BodyPool bodies;
    sqp::IslandPool islands, unrelated;
    const auto body = bodies.add_rigid();
    const auto island = islands.create();
    unrelated.create();
    REQUIRE(islands.add_body(island, body, bodies));
    CHECK_FALSE(bodies.remove(body, unrelated));
    CHECK_FALSE(bodies.clear(unrelated));
    CHECK_FALSE(unrelated.add_body(unrelated.create(), body, bodies));
    CHECK(bodies.is_valid(body));
    CHECK(bodies.island_of(body) == island);
    CHECK(islands.has_body(island, body));
}
