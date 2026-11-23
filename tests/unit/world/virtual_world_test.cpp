#include <limits>
#include <vector>

#include "sqphysics2d/detail/core_world.h"
#include <catch_amalgamated.hpp>

namespace {

class TestWorld : public sqp::VirtualWorld {
public:
    using sqp::VirtualWorld::install_stage;
    using sqp::VirtualWorld::set_entry_stage;

    std::vector<int> trace;

    StageId left_id = Finish;
    StageId right_id = Finish;
    StageId commit_id = Finish;
    StageId loop_id = Finish;
    StageId next_id = Finish;

    bool choose_right = false;
    int iterations = 0;

    static StageId start(VirtualWorld& base) {
        auto& world = static_cast<TestWorld&>(base);
        world.trace.push_back(1);
        return world.choose_right ? world.right_id : world.left_id;
    }

    static StageId left(VirtualWorld& base) {
        auto& world = static_cast<TestWorld&>(base);
        world.trace.push_back(2);
        return world.commit_id;
    }

    static StageId right(VirtualWorld& base) {
        auto& world = static_cast<TestWorld&>(base);
        world.trace.push_back(3);
        return world.commit_id;
    }

    static StageId commit(VirtualWorld& base) {
        auto& world = static_cast<TestWorld&>(base);
        world.trace.push_back(4);
        return Finish;
    }

    static StageId jump(VirtualWorld& base) {
        auto& world = static_cast<TestWorld&>(base);
        world.trace.push_back(5);
        return world.next_id;
    }

    static StageId loop(VirtualWorld& base) {
        auto& world = static_cast<TestWorld&>(base);
        ++world.iterations;
        return world.iterations < 4 ? world.loop_id : Finish;
    }
};
} // namespace

TEST_CASE("VirtualWorld stage registration", "[virtual_world]") {
    TestWorld world;

    const auto a = world.install_stage(TestWorld::start, "Start");
    const auto b = world.install_stage(TestWorld::left, "Left");
    const auto c = world.install_stage(TestWorld::commit, "Commit");

    REQUIRE(a == 2);
    REQUIRE(b == 3);
    REQUIRE(c == 4);
}

TEST_CASE("VirtualWorld state machine execution", "[virtual_world]") {
    TestWorld world;

    const auto start = world.install_stage(TestWorld::start, "Start");
    world.left_id = world.install_stage(TestWorld::left, "Left");
    world.right_id = world.install_stage(TestWorld::right, "Right");
    world.commit_id = world.install_stage(TestWorld::commit, "Commit");

    world.set_entry_stage(start);

    SECTION("Left branch") {
        REQUIRE(world.step());
        REQUIRE(world.trace == std::vector<int>({1, 2, 4}));
    }

    SECTION("Right branch") {
        world.choose_right = true;

        REQUIRE(world.step());
        REQUIRE(world.trace == std::vector<int>({1, 3, 4}));
    }

    SECTION("Repeated step") {
        REQUIRE(world.step());
        REQUIRE(world.step());

        REQUIRE(world.trace == std::vector<int>({1, 2, 4, 1, 2, 4}));
    }

    SECTION("Last registered stage is executable") {
        world.set_entry_stage(world.commit_id);

        REQUIRE(world.step());
        REQUIRE(world.trace == std::vector<int>({4}));
    }
}

TEST_CASE("VirtualWorld terminal stages", "[virtual_world]") {
    TestWorld world;

    SECTION("Default entry is Finish") {
        REQUIRE(world.step());
        REQUIRE(world.trace.empty());
    }

    SECTION("Explicit Finish") {
        world.set_entry_stage(sqp::VirtualWorld::Finish);
        REQUIRE(world.step());
    }

    SECTION("Explicit Error") {
        world.set_entry_stage(sqp::VirtualWorld::Error);
        REQUIRE_FALSE(world.step());
        CHECK(world.last_error_message() == "dispatcher: stage returned Error");
    }

    SECTION("Stage returns Error") {
        const auto id = world.install_stage(TestWorld::jump, "Jump");
        world.next_id = sqp::VirtualWorld::Error;
        world.set_entry_stage(id);

        REQUIRE_FALSE(world.step());
        REQUIRE(world.trace == std::vector<int>({5}));
        CHECK(world.last_error_message() == "Jump: stage returned Error");
    }

    SECTION("Stage returns Finish") {
        const auto id = world.install_stage(TestWorld::jump, "Jump");
        world.next_id = sqp::VirtualWorld::Finish;
        world.set_entry_stage(id);

        REQUIRE(world.step());
        REQUIRE(world.trace == std::vector<int>({5}));
    }
}

TEST_CASE("VirtualWorld invalid stage IDs", "[virtual_world]") {
    TestWorld world;
    const auto id = world.install_stage(TestWorld::jump, "Jump");

    SECTION("Invalid entry stage") {
        world.set_entry_stage(id + 1);
        REQUIRE_FALSE(world.step());
        CHECK(world.last_error_message() == "dispatcher: invalid stage ID 3");
    }

    SECTION("Invalid next stage") {
        world.next_id = id + 1;
        world.set_entry_stage(id);

        REQUIRE_FALSE(world.step());
        REQUIRE(world.trace == std::vector<int>({5}));
        CHECK(world.last_error_message() == "Jump: invalid stage ID 3");
    }

    SECTION("Large invalid stage ID") {
        world.next_id = std::numeric_limits<sqp::VirtualWorld::StageId>::max() - 2;
        world.set_entry_stage(id);

        REQUIRE_FALSE(world.step());
    }
}

TEST_CASE("VirtualWorld stage jumps and loops", "[virtual_world]") {
    TestWorld world;

    const auto jump = world.install_stage(TestWorld::jump, "Jump");
    world.loop_id = world.install_stage(TestWorld::loop, "Loop");

    SECTION("Forward jump") {
        world.next_id = world.loop_id;
        world.set_entry_stage(jump);

        REQUIRE(world.step());
        REQUIRE(world.iterations == 4);
    }

    SECTION("Self loop with termination") {
        world.set_entry_stage(world.loop_id);

        REQUIRE(world.step());
        REQUIRE(world.iterations == 4);
    }

    SECTION("Backward jump") {
        world.next_id = world.loop_id;
        world.set_entry_stage(world.loop_id);

        REQUIRE(world.step());
        REQUIRE(world.iterations == 4);
    }
}
