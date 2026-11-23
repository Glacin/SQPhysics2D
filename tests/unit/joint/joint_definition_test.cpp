#include <limits>
#include <type_traits>

#include "sqphysics2d/detail/data_storage.h"
#include <catch_amalgamated.hpp>

namespace {
template <class Def> void check_common_validation(sqp::BodyPool& bodies) {
    const auto a = bodies.add_rigid();
    const auto b = bodies.add_rigid();
    Def def;
    def.body_a = a;
    def.body_b = b;
    CHECK(sqp::validate_joint_def(def, bodies));

    def.body_a = {};
    CHECK(sqp::validate_joint_def(def, bodies));
    def.local_anchor_a = {1, 0};
    CHECK_FALSE(sqp::validate_joint_def(def, bodies));
    def.local_anchor_a = {};
    def.body_a = b;
    CHECK_FALSE(sqp::validate_joint_def(def, bodies));
    def.body_a = a;
    def.body_b = {};
    CHECK_FALSE(sqp::validate_joint_def(def, bodies));
    def.body_b = b;
    def.local_anchor_b.x = std::numeric_limits<sqp::real>::infinity();
    CHECK_FALSE(sqp::validate_joint_def(def, bodies));
    def.local_anchor_b = {};
    def.local_anchor_a.y = std::numeric_limits<sqp::real>::quiet_NaN();
    CHECK_FALSE(sqp::validate_joint_def(def, bodies));
    def.local_anchor_a = {};

    REQUIRE(bodies.remove(a));
    const auto replacement = bodies.add_rigid();
    CHECK(replacement.index() == a.index());
    CHECK_FALSE(sqp::validate_joint_def(def, bodies));
}
} // namespace

TEST_CASE("Joint definitions and runtime records have prescribed defaults", "[joint][definition]") {
    static_assert(std::is_base_of_v<sqp::JointDef, sqp::PinJointDef>);
    static_assert(std::is_base_of_v<sqp::JointDef, sqp::SliderJointDef>);
    static_assert(std::is_base_of_v<sqp::JointDef, sqp::WeldJointDef>);
    static_assert(std::is_base_of_v<sqp::JointDef, sqp::DistanceJointDef>);
    static_assert(std::is_base_of_v<sqp::JointDef, sqp::SpringJointDef>);
    CHECK(sqp::JointDef{}.body_a == sqp::BodyHandle{});
    CHECK(sqp::JointDef{}.body_b == sqp::BodyHandle{});
    CHECK_FALSE(sqp::JointDef{}.collide_connected);
    CHECK(sqp::SliderJointDef{}.local_axis_a == sqp::Vec2{1, 0});
    CHECK_FALSE(sqp::SliderJointDef{}.enable_limit);
    CHECK(sqp::SpringJointDef{}.hertz == 5);
    CHECK(sqp::SpringJointDef{}.damping_ratio == sqp::real(0.7));
    CHECK(sqp::SliderJoint{}.reference_angle == 0);
    CHECK(sqp::WeldJoint{}.reference_angle == 0);
    CHECK(sqp::DistanceJoint{}.length == 0);
    CHECK(sqp::SpringJoint{}.rest_length == 0);
    CHECK(sqp::JointRegistration{}.edges[0].prev == std::numeric_limits<sqp::u32>::max());
    CHECK(sqp::JointRegistration{}.edges[1].next == std::numeric_limits<sqp::u32>::max());
}

TEST_CASE("Joint handles use sortable slots and generations", "[joint][handle]") {
    constexpr sqp::JointHandle empty{};
    constexpr sqp::JointHandle first{0, 1};
    constexpr sqp::JointHandle recycled{0, 2};
    constexpr sqp::JointHandle next{1, 1};
    static_assert(!empty.valid());
    static_assert(first.valid() && first < recycled && recycled < next);
    static_assert(!sqp::JointHandle{sqp::max_joint_slot + 1, 1}.valid());
    static_assert(!sqp::JointHandle{0, 0}.valid());
    static_assert(sqp::joint_edge_key(first, 0) == 0);
    static_assert(sqp::joint_edge_key(first, 1) == 1);
    static_assert(sqp::joint_edge_key(next, 0) == 2);
    static_assert(sqp::joint_edge_key({sqp::max_joint_slot, 1}, 1) == 0x7fffffffu);
    static_assert(sqp::joint_edge_key(empty, 0) == std::numeric_limits<sqp::u32>::max());
    static_assert(sqp::joint_edge_key(first, 2) == std::numeric_limits<sqp::u32>::max());
    CHECK(first != recycled);
}

TEST_CASE("Joint definitions reject invalid endpoints", "[joint][definition]") {
    sqp::BodyPool bodies;
    check_common_validation<sqp::PinJointDef>(bodies);
    check_common_validation<sqp::SliderJointDef>(bodies);
    check_common_validation<sqp::WeldJointDef>(bodies);
    check_common_validation<sqp::DistanceJointDef>(bodies);
    check_common_validation<sqp::SpringJointDef>(bodies);
}

TEST_CASE("Slider validates its axis and limits even when disabled", "[joint][definition]") {
    sqp::BodyPool bodies;
    sqp::SliderJointDef def;
    def.body_b = bodies.add_rigid();
    def.local_axis_a = {3, 4};
    CHECK(sqp::validate_joint_def(def, bodies));
    CHECK(def.local_axis_a == sqp::Vec2{3, 4});
    def.local_axis_a = {};
    CHECK_FALSE(sqp::validate_joint_def(def, bodies));
    def.local_axis_a = {1, 0};
    def.lower_translation = 2;
    def.upper_translation = 1;
    CHECK_FALSE(sqp::validate_joint_def(def, bodies));
    def.upper_translation = 2;
    CHECK(sqp::validate_joint_def(def, bodies));
    def.upper_translation = std::numeric_limits<sqp::real>::infinity();
    CHECK_FALSE(sqp::validate_joint_def(def, bodies));
}

TEST_CASE("Distance and Spring validate scalar parameters", "[joint][definition]") {
    sqp::BodyPool bodies;
    const auto b = bodies.add_rigid();
    sqp::DistanceJointDef distance;
    distance.body_b = b;
    CHECK(sqp::validate_joint_def(distance, bodies));
    distance.length = -1;
    CHECK_FALSE(sqp::validate_joint_def(distance, bodies));
    distance.length = std::numeric_limits<sqp::real>::quiet_NaN();
    CHECK_FALSE(sqp::validate_joint_def(distance, bodies));

    sqp::SpringJointDef spring;
    spring.body_b = b;
    spring.rest_length = 0;
    CHECK(sqp::validate_joint_def(spring, bodies));
    spring.rest_length = -1;
    CHECK_FALSE(sqp::validate_joint_def(spring, bodies));
    spring.rest_length = 0;
    spring.hertz = 0;
    CHECK_FALSE(sqp::validate_joint_def(spring, bodies));
    spring.hertz = 5;
    spring.damping_ratio = -1;
    CHECK_FALSE(sqp::validate_joint_def(spring, bodies));
    spring.damping_ratio = std::numeric_limits<sqp::real>::infinity();
    CHECK_FALSE(sqp::validate_joint_def(spring, bodies));
}
