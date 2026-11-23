#include <algorithm>
#include <random>
#include <vector>

#include "sqphysics2d/detail/broad_phase.h"
#include <catch_amalgamated.hpp>

TEST_CASE("dynamic Aabb tree is usable independently of broad_phase policy",
          "[broad_phase][dynamic-tree]") {

    sqp::DynamicAabbTree tree;

    tree.add_proxy(5, {{0, 0}, {2, 2}}, sqp::BroadPhase::PROXY_DYNAMIC);
    tree.add_proxy(1, {{1, 1}, {3, 3}}, sqp::BroadPhase::PROXY_STATIC);

    std::vector<sqp::BroadPhase::Pair> pairs;
    tree.build_pairs(pairs);

    REQUIRE(pairs.size() == 1);
    CHECK(pairs[0].a == 1);
    CHECK(pairs[0].b == 5);
}

TEST_CASE("three scalar trees preserve SQP cross-type pairs through local rebuild",
          "[broad_phase][rebuild]") {
    sqp::DynamicAabbTree tree;
    const auto fixed = tree.add_proxy(0, {{0, 0}, {2, 2}}, sqp::BroadPhase::PROXY_STATIC);
    const auto moving = tree.add_proxy(1, {{1, 0}, {3, 2}}, sqp::BroadPhase::PROXY_DYNAMIC);
    const auto driven = tree.add_proxy(2, {{1, 1}, {3, 3}}, sqp::BroadPhase::PROXY_KINEMATIC);
    const auto driven_peer =
        tree.add_proxy(3, {{1, 1}, {3, 3}}, sqp::BroadPhase::PROXY_KINEMATIC);
    tree.add_proxy(4, {{0, 0}, {2, 2}}, sqp::BroadPhase::PROXY_STATIC);

    std::vector<sqp::BroadPhase::Pair> before, after;
    tree.build_pairs(before);
    CHECK(before.size() == 9);
    CHECK(std::find(before.begin(), before.end(), sqp::BroadPhase::Pair{0, 2}) !=
          before.end());
    CHECK(std::find(before.begin(), before.end(), sqp::BroadPhase::Pair{2, 3}) !=
          before.end());
    tree.rebuild_changed();
    tree.build_pairs(after);
    CHECK(after == before);
    CHECK(tree.proxy_matches(fixed, 0));
    CHECK(tree.proxy_matches(moving, 1));
    CHECK(tree.proxy_matches(driven, 2));
    CHECK(tree.proxy_matches(driven_peer, 3));

    REQUIRE(tree.update_proxy(driven, {{10, 10}, {11, 11}},
                              sqp::BroadPhase::PROXY_DYNAMIC));
    tree.rebuild_changed();
    CHECK(tree.proxy_matches(driven, 2));
    std::vector<sqp::u32> hits;
    tree.query_aabb({{10, 10}, {11, 11}}, hits);
    CHECK(hits == std::vector<sqp::u32>{2});
}

TEST_CASE("Aabb helpers merge and sweep bounds", "[broad_phase][aabb]") {

    const sqp::Aabb a{{-2, -1}, {1, 3}};
    const sqp::Aabb b{{0, -4}, {5, 2}};
    const sqp::Aabb merged = sqp::Aabb::merge(a, b);

    CHECK(merged.min.x == Catch::Approx(-2));
    CHECK(merged.min.y == Catch::Approx(-4));
    CHECK(merged.max.x == Catch::Approx(5));
    CHECK(merged.max.y == Catch::Approx(3));

    const sqp::Aabb swept = sqp::Aabb::swept(a, {3, -2});
    CHECK(swept.min.y == Catch::Approx(-3));
    CHECK(swept.max.x == Catch::Approx(4));
}

TEST_CASE("broad_phase reports normalized overlaps and skips static pairs", "[broad_phase]") {

    sqp::DynamicAabbTree broad_phase;

    broad_phase.add_proxy(4, {{0, 0}, {2, 2}}, sqp::BroadPhase::PROXY_STATIC);
    broad_phase.add_proxy(2, {{1, 1}, {3, 3}}, sqp::BroadPhase::PROXY_STATIC);
    broad_phase.add_proxy(3, {{sqp::real(1.5), 0}, {sqp::real(2.5), 2}});

    std::vector<sqp::BroadPhase::Pair> pairs;
    broad_phase.build_pairs(pairs);

    REQUIRE(pairs.size() == 2);
    CHECK(pairs[0].a == 2);
    CHECK(pairs[0].b == 3);
    CHECK(pairs[1].a == 3);
    CHECK(pairs[1].b == 4);
}

TEST_CASE("broad_phase update removes the old tree path", "[broad_phase][update]") {

    sqp::DynamicAabbTree broad_phase;

    const sqp::i32 handle = broad_phase.add_proxy(1, {{0, 0}, {1, 1}});
    broad_phase.add_proxy(2, {{sqp::real(0.5), 0}, {sqp::real(1.5), 1}});
    broad_phase.add_proxy(3, {{10, 10}, {12, 12}});

    REQUIRE(broad_phase.update_proxy(handle, {{sqp::real(10.5), sqp::real(10.5)}, {11, 11}},
                                     sqp::BroadPhase::PROXY_DYNAMIC));
    std::vector<sqp::u32> ids;

    broad_phase.query_aabb({{0, 0}, {2, 2}}, ids);
    CHECK(ids == std::vector<sqp::u32>{2});

    broad_phase.query_aabb({{10, 10}, {12, 12}}, ids);
    CHECK(ids == std::vector<sqp::u32>{1, 3});
}

TEST_CASE("querying pairs by proxy handle uses stored bounds and flags", "[broad_phase][pair]") {
    sqp::DynamicAabbTree tree;
    const auto selected = tree.add_proxy(0, {{0, 0}, {2, 2}}, sqp::BroadPhase::PROXY_STATIC);
    const auto static_peer = tree.add_proxy(1, {{1, 1}, {3, 3}}, sqp::BroadPhase::PROXY_STATIC);
    tree.add_proxy(2, {{1, 1}, {3, 3}}, sqp::BroadPhase::PROXY_DYNAMIC);

    std::vector<sqp::BroadPhase::Pair> by_bounds;
    std::vector<sqp::BroadPhase::Pair> by_handle;
    tree.query_pairs(0, {{0, 0}, {2, 2}}, sqp::BroadPhase::PROXY_STATIC, by_bounds);
    REQUIRE(tree.query_proxy_pairs(selected, by_handle));
    CHECK(by_handle == by_bounds);
    CHECK(by_handle == std::vector<sqp::BroadPhase::Pair>{{0, 2}});

    REQUIRE(tree.update_proxy(selected, {{0, 0}, {2, 2}}, sqp::BroadPhase::PROXY_DYNAMIC));
    by_handle.clear();
    REQUIRE(tree.query_proxy_pairs(selected, by_handle));
    CHECK(by_handle.size() == 2);
    CHECK(std::find(by_handle.begin(), by_handle.end(), sqp::BroadPhase::Pair{0, 1}) !=
          by_handle.end());
    CHECK(std::find(by_handle.begin(), by_handle.end(), sqp::BroadPhase::Pair{0, 2}) !=
          by_handle.end());

    REQUIRE(tree.update_proxy(selected, {{8, 8}, {9, 9}}, sqp::BroadPhase::PROXY_DYNAMIC));
    by_bounds.clear();
    by_handle.clear();
    tree.query_pairs(0, {{8, 8}, {9, 9}}, sqp::BroadPhase::PROXY_DYNAMIC, by_bounds);
    REQUIRE(tree.query_proxy_pairs(selected, by_handle));
    CHECK(by_handle == by_bounds);
    CHECK(by_handle.empty());

    CHECK_FALSE(tree.query_proxy_pairs(-1, by_handle));
    REQUIRE(tree.remove_proxy(static_peer));
    CHECK_FALSE(tree.query_proxy_pairs(static_peer, by_handle));
    CHECK(by_handle.empty());
}

TEST_CASE("shrinking a broad_phase proxy removes its previous extent", "[broad_phase][update]") {

    sqp::DynamicAabbTree broad_phase;

    const sqp::i32 handle = broad_phase.add_proxy(9, {{0, 0}, {10, 10}});
    REQUIRE(broad_phase.update_proxy(handle, {{4, 4}, {6, 6}}, sqp::BroadPhase::PROXY_DYNAMIC));

    std::vector<sqp::u32> ids;

    broad_phase.query_aabb({{0, 0}, {1, 1}}, ids);
    CHECK(ids.empty());

    broad_phase.query_aabb({{5, 5}, {5, 5}}, ids);
    CHECK(ids == std::vector<sqp::u32>{9});
}

TEST_CASE("exact swept tree supports stable removal", "[broad_phase][exact][remove][oracle]") {

    sqp::DynamicAabbTree tree;
    const sqp::i32 first = tree.add_proxy(1, {{0, 0}, {1, 1}});
    const sqp::i32 second = tree.add_proxy(2, {{sqp::real(1.05), 0}, {2, 1}});

    bool reinserted = true;
    REQUIRE(tree.update_proxy(first, {{sqp::real(0.02), 0}, {sqp::real(1.02), 1}},
                              sqp::BroadPhase::PROXY_DYNAMIC, &reinserted));
    CHECK_FALSE(reinserted);

    std::vector<sqp::BroadPhase::Pair> pairs;
    tree.build_pairs(pairs);
    CHECK(pairs.empty());

    std::vector<sqp::u32> ids;
    tree.query_aabb({{sqp::real(1.03), 0}, {sqp::real(1.04), 1}}, ids);
    CHECK(ids.empty());

    REQUIRE(tree.remove_proxy(second));
    CHECK(tree.proxy_count() == 1);
    CHECK_FALSE(tree.remove_proxy(second));
    tree.query_aabb({{1, 0}, {2, 1}}, ids);
    CHECK(ids == std::vector<sqp::u32>{1});
}

TEST_CASE("fat tree bounds retain exact pair and ray results", "[broad_phase][fat][oracle]") {
    sqp::DynamicAabbTree tree;
    const auto moving = tree.add_proxy(1, {{0, 0}, {1, 1}});
    const auto peer = tree.add_proxy(2, {{sqp::real(1.05), 0}, {2, 1}});
    bool reinserted = true;
    REQUIRE(tree.update_proxy(moving, {{sqp::real(0.02), 0}, {sqp::real(1.02), 1}},
                              sqp::BroadPhase::PROXY_DYNAMIC, &reinserted));
    CHECK_FALSE(reinserted);

    std::vector<sqp::BroadPhase::Pair> pairs;
    REQUIRE(tree.query_proxy_pairs(moving, pairs));
    CHECK(pairs.empty());
    tree.query_pairs(99, {{sqp::real(1.03), 0}, {sqp::real(1.04), 1}},
                     sqp::BroadPhase::PROXY_DYNAMIC, pairs);
    CHECK(pairs.empty());
    tree.build_pairs(pairs);
    CHECK(pairs.empty());

    std::vector<sqp::u32> ids;
    tree.query_aabb({{sqp::real(1.03), 0}, {sqp::real(1.04), 1}}, ids);
    CHECK(ids.empty());

    struct RayContext { std::vector<sqp::u32> ids; } context;
    const auto callback = [](const sqp::BroadPhase::RayCastInput&, sqp::u32 id,
                             void* opaque) -> sqp::real {
        static_cast<RayContext*>(opaque)->ids.push_back(id);
        return sqp::real(-1);
    };
    tree.ray_cast({{sqp::real(1.03), -1}, {0, 3}, 1}, callback, &context);
    CHECK(context.ids.empty());
    tree.ray_cast({{sqp::real(0.95), -1}, {0, 3}, 1}, callback, &context);
    CHECK(context.ids == std::vector<sqp::u32>{1});

    REQUIRE(tree.update_proxy(moving, {{sqp::real(1.5), 0}, {sqp::real(2.5), 1}},
                              sqp::BroadPhase::PROXY_DYNAMIC, &reinserted));
    CHECK(reinserted);
    tree.build_pairs(pairs);
    CHECK(pairs == std::vector<sqp::BroadPhase::Pair>{{1, 2}});

    REQUIRE(tree.remove_proxy(peer));
    tree.rebuild();
    tree.build_pairs(pairs);
    CHECK(pairs.empty());
}

TEST_CASE("dynamic proxy displacement amortizes reinsertion while static bounds stay exact",
          "[broad_phase][fat][update]") {
    sqp::DynamicAabbTree tree;
    tree.set_proxy_extension_ratio(sqp::real(0.5));
    const auto moving = tree.add_proxy(1, {{0, 0}, {1, 1}});
    const auto fixed = tree.add_proxy(2, {{3, 0}, {4, 1}}, sqp::BroadPhase::PROXY_STATIC);
    unsigned reinsertion_count = 0;
    for (int step = 1; step <= 20; ++step) {
        const sqp::real offset = sqp::real(step) * sqp::real(0.078125);
        bool reinserted = false;
        REQUIRE(tree.update_proxy(moving, {{0, offset}, {1, 1 + offset}},
                                  sqp::BroadPhase::PROXY_DYNAMIC, &reinserted));
        reinsertion_count += reinserted ? 1u : 0u;
    }
    CHECK(reinsertion_count < 10);

    bool reinserted = false;
    REQUIRE(tree.update_proxy(fixed, {{sqp::real(3.1), 0}, {sqp::real(4.1), 1}},
                              sqp::BroadPhase::PROXY_STATIC, &reinserted));
    CHECK(reinserted);
    std::vector<sqp::u32> ids;
    tree.query_aabb({{3, 0}, {sqp::real(3.05), 1}}, ids);
    CHECK(ids.empty());
}

TEST_CASE("proxy extension uses each swept axis and live ratio changes keep handles",
          "[broad_phase][fat][settings]") {
    sqp::DynamicAabbTree tree;
    tree.set_proxy_extension_ratio(sqp::real(0.2));
    const auto handle = tree.add_proxy(7, {{0, 0}, {2, 4}});
    bool reinserted = true;
    REQUIRE(tree.update_proxy(handle, {{sqp::real(0.19), 0}, {sqp::real(2.19), 4}},
                              sqp::BroadPhase::PROXY_DYNAMIC, &reinserted));
    CHECK_FALSE(reinserted);
    REQUIRE(tree.update_proxy(handle, {{sqp::real(0.19), sqp::real(0.41)},
                                       {sqp::real(2.19), sqp::real(4.41)}},
                              sqp::BroadPhase::PROXY_DYNAMIC, &reinserted));
    CHECK(reinserted);

    tree.set_proxy_extension_ratio(0);
    CHECK(tree.proxy_matches(handle, 7));
    std::vector<sqp::u32> ids;
    tree.query_aabb({{sqp::real(0.19), sqp::real(0.41)},
                     {sqp::real(2.19), sqp::real(4.41)}}, ids);
    CHECK(ids == std::vector<sqp::u32>{7});
    REQUIRE(tree.update_proxy(handle, {{sqp::real(0.2), sqp::real(0.41)},
                                       {sqp::real(2.2), sqp::real(4.41)}},
                              sqp::BroadPhase::PROXY_DYNAMIC, &reinserted));
    CHECK(reinserted);
}

TEST_CASE("tree updates agree with a brute force overlap oracle", "[broad_phase][oracle]") {

    constexpr sqp::u32 count = 32;
    std::mt19937 rng(0x05162020);
    std::uniform_real_distribution<float> position(-20.0f, 20.0f);

    std::vector<sqp::Aabb> bounds(count);
    std::vector<sqp::i32> handles(count);

    sqp::DynamicAabbTree broad_phase;
    for (sqp::u32 i = 0; i < count; ++i) {
        const sqp::real x = static_cast<sqp::real>(position(rng));
        const sqp::real y = static_cast<sqp::real>(position(rng));
        bounds[i] = {{x, y}, {x + sqp::real(2), y + sqp::real(2)}};
        handles[i] = broad_phase.add_proxy(i, bounds[i]);
    }

    bool all_updates_succeeded = true;
    bool all_queries_matched = true;

    sqp::u32 first_failed_step = 0;

    for (sqp::u32 step = 0; step < 64; ++step) {
        const sqp::u32 moved = step % count;
        const sqp::real x = static_cast<sqp::real>(position(rng));
        const sqp::real y = static_cast<sqp::real>(position(rng));

        bounds[moved] = {{x, y}, {x + sqp::real(2), y + sqp::real(2)}};

        const bool updated =
            broad_phase.update_proxy(handles[moved], bounds[moved], sqp::BroadPhase::PROXY_DYNAMIC);

        if (!updated && all_updates_succeeded)
            first_failed_step = step;

        all_updates_succeeded = updated && all_updates_succeeded;
        const sqp::Aabb query{{x - sqp::real(1), y - sqp::real(1)},
                              {x + sqp::real(3), y + sqp::real(3)}};
        std::vector<sqp::u32> expected;

        for (sqp::u32 i = 0; i < count; ++i) {
            if (sqp::Aabb::overlap(query, bounds[i])) {
                expected.push_back(i);
            }
        }

        std::vector<sqp::u32> actual;

        broad_phase.query_aabb(query, actual);
        if (actual != expected && all_queries_matched)
            first_failed_step = step;

        all_queries_matched = actual == expected && all_queries_matched;
        std::array<sqp::DynamicAabbTree::QualityStats, 3> quality;
        CHECK(broad_phase.inspect(quality));
        CHECK(quality[2].leaf_count == count);
    }

    CAPTURE(first_failed_step);
    REQUIRE(all_updates_succeeded);
    CHECK(all_queries_matched);
}

TEST_CASE("broad_phase bulk rebuild preserves handles and overlap results",
          "[broad_phase][rebuild][oracle]") {

    constexpr sqp::u32 count = 127;

    sqp::DynamicAabbTree tree;
    std::vector<sqp::Aabb> bounds(count);
    std::vector<sqp::i32> handles(count);

    for (sqp::u32 i = 0; i < count; ++i) {
        const sqp::real x = static_cast<sqp::real>((i * 37) % 53);
        const sqp::real y = static_cast<sqp::real>((i * 19) % 47);

        bounds[i] = {{x, y}, {x + sqp::real(3), y + sqp::real(2)}};
        handles[i] = tree.add_proxy(i, bounds[i]);
    }

    std::vector<sqp::BroadPhase::Pair> before;
    std::vector<sqp::BroadPhase::Pair> after;

    tree.build_pairs(before);
    tree.rebuild();
    tree.build_pairs(after);
    CHECK(before == after);

    bounds[63] = {{-10, -10}, {-8, -8}};
    REQUIRE(tree.update_proxy(handles[63], bounds[63], sqp::BroadPhase::PROXY_DYNAMIC));
    std::vector<sqp::u32> ids;
    tree.query_aabb({{-11, -11}, {-7, -7}}, ids);
    CHECK(ids == std::vector<sqp::u32>{63});

    tree.rebuild();
    tree.query_aabb({{-11, -11}, {-7, -7}}, ids);
    CHECK(ids == std::vector<sqp::u32>{63});
}

TEST_CASE("partial SAH rebuild preserves pairs through moves and removals",
          "[broad_phase][rebuild][oracle]") {
    sqp::DynamicAabbTree tree;
    struct Entry {
        sqp::i32 handle;
        sqp::u32 id;
        sqp::Aabb bounds;
    };
    std::vector<Entry> entries;
    for (sqp::u32 i = 0; i < 40; ++i) {
        const sqp::real x = static_cast<sqp::real>(i % 7);
        const sqp::real y = static_cast<sqp::real>(i % 5);
        const sqp::Aabb box{{x, y}, {x + 2, y + 2}};
        entries.push_back({tree.add_proxy(i, box), i, box});
    }
    for (sqp::u32 step = 0; step < 80; ++step) {
        Entry& changed = entries[step % entries.size()];
        if (step % 9 == 0) {
            REQUIRE(tree.remove_proxy(changed.handle));
            changed.id = 40 + step;
            changed.handle = tree.add_proxy(changed.id, changed.bounds);
        } else {
            const sqp::real x = static_cast<sqp::real>((step * 11) % 17);
            const sqp::real y = static_cast<sqp::real>((step * 7) % 13);
            changed.bounds = {{x, y}, {x + 2, y + 2}};
            REQUIRE(tree.update_proxy(changed.handle, changed.bounds,
                                      sqp::BroadPhase::PROXY_DYNAMIC));
        }
        tree.rebuild_changed();
        std::array<sqp::DynamicAabbTree::QualityStats, 3> quality;
        REQUIRE(tree.inspect(quality));
        CHECK(quality[2].leaf_count == entries.size());
        std::vector<sqp::BroadPhase::Pair> expected, actual;
        for (std::size_t a = 0; a < entries.size(); ++a)
            for (std::size_t b = a + 1; b < entries.size(); ++b)
                if (sqp::Aabb::overlap(entries[a].bounds, entries[b].bounds))
                    expected.push_back({std::min(entries[a].id, entries[b].id),
                                        std::max(entries[a].id, entries[b].id)});
        std::sort(expected.begin(), expected.end());
        tree.build_pairs(actual);
        CHECK(actual == expected);
        for (const Entry& entry : entries)
            CHECK(tree.proxy_matches(entry.handle, entry.id));
    }
}

TEST_CASE("broad_phase pair traversal handles coincident leaves deterministically",
          "[broad_phase][rebuild]") {

    sqp::DynamicAabbTree tree;

    for (sqp::u32 i = 0; i < 32; ++i) {
        tree.add_proxy(31 - i, {{0, 0}, {1, 1}});
    }
    tree.rebuild();
    std::array<sqp::DynamicAabbTree::QualityStats, 3> quality;
    REQUIRE(tree.inspect(quality));
    CHECK(quality[2].leaf_count == 32);

    std::vector<sqp::BroadPhase::Pair> first;
    std::vector<sqp::BroadPhase::Pair> second;

    tree.build_pairs(first);
    tree.rebuild();
    tree.build_pairs(second);

    CHECK(first.size() == 32u * 31u / 2u);
    CHECK(first == second);
}

TEST_CASE("unbuilt broad_phase batch becomes queryable after rebuild", "[broad_phase][rebuild]") {

    sqp::DynamicAabbTree tree;

    tree.add_proxy_unbuilt(7, {{0, 0}, {2, 2}});
    tree.add_proxy_unbuilt(3, {{1, 1}, {3, 3}});

    std::vector<sqp::BroadPhase::Pair> pairs;
    tree.build_pairs(pairs);

    CHECK(pairs.empty());
    tree.rebuild();
    tree.build_pairs(pairs);

    REQUIRE(pairs.size() == 1);
    CHECK(pairs[0].a == 3);
    CHECK(pairs[0].b == 7);
}

TEST_CASE("dynamic tree raycast supports clipping skipping and termination",
          "[broad_phase][raycast]") {

    sqp::DynamicAabbTree tree;
    tree.add_proxy(1, {{0, -1}, {1, 1}});
    tree.add_proxy(2, {{4, -1}, {5, 1}});
    tree.add_proxy(3, {{8, -1}, {9, 1}});

    struct Context {
        std::vector<sqp::u32> ids;
        bool terminate = false;
    } context;

    const auto callback = [](const sqp::BroadPhase::RayCastInput&, sqp::u32 id,
                             void* opaque) -> sqp::real {
        auto& context = *static_cast<Context*>(opaque);
        context.ids.push_back(id);
        if (context.terminate)
            return sqp::real(0);
        if (id == 1)
            return sqp::real(0.25);
        return sqp::real(-1);
    };

    std::vector<sqp::i32> stack;
    const auto stats = tree.ray_cast({{-2, 0}, {12, 0}, 1}, callback, &context, &stack);
    CHECK(context.ids == std::vector<sqp::u32>{1});
    CHECK(stats.nodes_visited > 0);
    CHECK(stats.leaves_visited == 1);

    context = {};
    context.terminate = true;
    const auto terminated = tree.ray_cast({{-2, 0}, {12, 0}, 1}, callback, &context, &stack);
    CHECK(context.ids.size() == 1);
    CHECK(terminated.leaves_visited == 1);

    context = {};
    CHECK(tree.ray_cast({{0, 0}, {}, 1}, callback, &context).nodes_visited == 0);
    CHECK(context.ids.empty());
}

TEST_CASE("broad_phase pair ordering is lexicographic", "[broad_phase][pair]") {

    CHECK(sqp::BroadPhase::Pair{1, 9} < sqp::BroadPhase::Pair{2, 0});
    CHECK(sqp::BroadPhase::Pair{2, 3} < sqp::BroadPhase::Pair{2, 4});
    CHECK_FALSE(sqp::BroadPhase::Pair{2, 4} < sqp::BroadPhase::Pair{2, 4});
}
