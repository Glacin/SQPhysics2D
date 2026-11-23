#include <vector>

#include "sqphysics2d/detail/data_storage.h"
#include <catch_amalgamated.hpp>

namespace {
std::vector<sqp::JointHandle> attached(const sqp::JointPool& joints, const sqp::BodyPool& bodies,
                                       sqp::BodyHandle body) {
    std::vector<sqp::JointHandle> result;
    for (auto joint = joints.first_attached(bodies, body); joint.valid();
         joint = joints.next_attached(body, joint)) {
        result.push_back(joint);
        REQUIRE(result.size() < 100);
    }
    return result;
}

template <class Def> void check_type_swap_remove() {
    sqp::BodyPool bodies;
    sqp::JointPool joints;
    Def def;
    def.body_a = bodies.add_rigid();
    def.body_b = bodies.add_rigid();
    const auto first = joints.create(bodies, def);
    const auto second = joints.create(bodies, def);
    REQUIRE(first.valid());
    REQUIRE(second.valid());
    REQUIRE(joints.destroy(bodies, first));
    REQUIRE(joints.get(second));
    CHECK(joints.get(second)->type_index == 0);
    CHECK(joints.is_valid(second));
    CHECK_FALSE(joints.is_valid(first));
}
} // namespace

TEST_CASE("Every joint type repairs its dense owner after swap removal", "[joint][pool]") {
    check_type_swap_remove<sqp::PinJointDef>();
    check_type_swap_remove<sqp::SliderJointDef>();
    check_type_swap_remove<sqp::WeldJointDef>();
    check_type_swap_remove<sqp::DistanceJointDef>();
    check_type_swap_remove<sqp::SpringJointDef>();
}

TEST_CASE("JointPool stores five kinds and keeps dense owners after removal", "[joint][pool]") {
    sqp::BodyPool bodies;
    sqp::JointPool joints;
    auto a = bodies.add_rigid();
    auto b = bodies.add_rigid();
    REQUIRE(bodies.set_pose(a, {}, 0.25f));
    REQUIRE(bodies.set_pose(b, {}, 1.5f));

    sqp::PinJointDef pin;
    pin.body_a = a;
    pin.body_b = b;
    pin.local_anchor_a = {1, 2};
    pin.collide_connected = true;
    sqp::SliderJointDef slider;
    slider.body_a = a;
    slider.body_b = b;
    slider.local_axis_a = {3, 4};
    slider.enable_limit = true;
    slider.lower_translation = -2;
    slider.upper_translation = 5;
    sqp::WeldJointDef weld;
    weld.body_a = a;
    weld.body_b = b;
    sqp::DistanceJointDef distance;
    distance.body_a = a;
    distance.body_b = b;
    distance.length = 7;
    sqp::SpringJointDef spring;
    spring.body_a = a;
    spring.body_b = b;
    spring.rest_length = 9;
    spring.hertz = 3;
    spring.damping_ratio = 0.4f;

    const auto p = joints.create(bodies, pin);
    const auto s = joints.create(bodies, slider);
    const auto w = joints.create(bodies, weld);
    const auto d = joints.create(bodies, distance);
    const auto sp = joints.create(bodies, spring);
    REQUIRE(p.valid());
    REQUIRE(s.valid());
    REQUIRE(w.valid());
    REQUIRE(d.valid());
    REQUIRE(sp.valid());
    REQUIRE(joints.get_pin(p));
    CHECK(joints.get_pin(p)->local_anchor_a == sqp::Vec2{1, 2});
    CHECK(joints.get(p)->collide_connected);
    REQUIRE(joints.get_slider(s));
    CHECK(joints.get_slider(s)->local_axis_a == sqp::Vec2{3, 4});
    CHECK(joints.get_slider(s)->reference_angle == Catch::Approx(1.25f));
    REQUIRE(joints.get_weld(w));
    CHECK(joints.get_weld(w)->reference_angle == Catch::Approx(1.25f));
    REQUIRE(joints.get_distance(d));
    CHECK(joints.get_distance(d)->length == 7);
    REQUIRE(joints.get_spring(sp));
    CHECK(joints.get_spring(sp)->rest_length == 9);
    CHECK(joints.get_pin(s) == nullptr);
    CHECK(attached(joints, bodies, a) == std::vector<sqp::JointHandle>{sp, d, w, s, p});
    CHECK(attached(joints, bodies, b) == std::vector<sqp::JointHandle>{sp, d, w, s, p});

    const auto p2 = joints.create(bodies, pin);
    REQUIRE(p2.valid());
    REQUIRE(joints.destroy(bodies, w));
    CHECK(attached(joints, bodies, a) == std::vector<sqp::JointHandle>{p2, sp, d, s, p});
    CHECK(joints.get(d)->edges[0].next == sqp::joint_edge_key(s, 0));
    CHECK(joints.get(s)->edges[0].prev == sqp::joint_edge_key(d, 0));
    REQUIRE(joints.destroy(bodies, p));
    REQUIRE(joints.get_pin(p2));
    CHECK(joints.get(p2)->type_index == 0);
    CHECK_FALSE(joints.is_valid(p));
    CHECK_FALSE(joints.destroy(bodies, p));
    const auto recycled = joints.create(bodies, pin);
    REQUIRE(recycled.valid());
    CHECK(recycled.index == p.index);
    CHECK(recycled.generation != p.generation);
    CHECK(joints.get_pin(p) == nullptr);
    CHECK_FALSE(bodies.remove(a));
}

TEST_CASE("JointPool world anchor and failed creation preserve adjacency", "[joint][pool]") {
    sqp::BodyPool bodies;
    sqp::JointPool joints;
    auto b = bodies.add_rigid();
    sqp::WeldJointDef world;
    world.body_b = b;
    const auto one = joints.create(bodies, world);
    const auto two = joints.create(bodies, world);
    REQUIRE(one.valid());
    REQUIRE(two.valid());
    CHECK(joints.get_weld(one)->reference_angle == 0);
    CHECK(joints.get(one)->edges[0].next == sqp::INVALID_JOINT_EDGE);
    CHECK(joints.get(one)->edges[0].prev == sqp::INVALID_JOINT_EDGE);
    CHECK(attached(joints, bodies, b) == std::vector<sqp::JointHandle>{two, one});
    CHECK(joints.get(two)->edges[1].next == sqp::joint_edge_key(one, 1));
    CHECK(joints.get(one)->edges[1].prev == sqp::joint_edge_key(two, 1));

    world.local_anchor_a = {1, 0};
    CHECK_FALSE(joints.create(bodies, world).valid());
    world.local_anchor_a = {};
    world.body_b = {};
    CHECK_FALSE(joints.create(bodies, world).valid());
    CHECK(attached(joints, bodies, b) == std::vector<sqp::JointHandle>{two, one});
    CHECK_FALSE(bodies.clear());
    REQUIRE(joints.destroy_attached(bodies, b));
    CHECK(attached(joints, bodies, b).empty());
    CHECK_FALSE(joints.is_valid(one));
    CHECK_FALSE(joints.is_valid(two));
    CHECK(bodies.clear());
}

TEST_CASE("Body dense moves preserve joint chains and cascading removal", "[joint][body]") {
    sqp::ShapePool shapes;
    sqp::CircleShapeDef circle;
    circle.radius = 1;
    sqp::PolygonShapeDef box;
    box.vertices = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
    sqp::BodyPool bodies;
    sqp::JointPool joints;
    const auto r0 = bodies.add_rigid();
    const auto r1 = bodies.add_rigid();
    const auto r2 = bodies.add_rigid();
    const auto c0 = bodies.add_lite_circle({}, shapes.add_circle(circle), shapes);
    const auto c1 = bodies.add_lite_circle({}, shapes.add_circle(circle), shapes);
    const auto x0 = bodies.add_aabb_box({}, shapes.add_polygon(box), shapes);
    const auto x1 = bodies.add_aabb_box({}, shapes.add_polygon(box), shapes);
    REQUIRE(bodies.is_valid(c0));
    REQUIRE(bodies.is_valid(c1));
    REQUIRE(bodies.is_valid(x0));
    REQUIRE(bodies.is_valid(x1));
    sqp::PinJointDef def;
    def.body_a = r2;
    def.body_b = c1;
    const auto rc = joints.create(bodies, def);
    def.body_a = x1;
    const auto xc = joints.create(bodies, def);
    def.body_a = r2;
    def.body_b = x1;
    const auto rx = joints.create(bodies, def);
    REQUIRE(rc.valid());
    REQUIRE(xc.valid());
    REQUIRE(rx.valid());

    REQUIRE(bodies.remove(r1));
    REQUIRE(bodies.remove(c0));
    REQUIRE(bodies.remove(x0));
    CHECK(attached(joints, bodies, r2) == std::vector<sqp::JointHandle>{rx, rc});
    CHECK(attached(joints, bodies, c1) == std::vector<sqp::JointHandle>{xc, rc});
    CHECK(attached(joints, bodies, x1) == std::vector<sqp::JointHandle>{rx, xc});
    CHECK_FALSE(bodies.remove(c1));
    REQUIRE(bodies.remove(c1, joints));
    CHECK_FALSE(joints.is_valid(rc));
    CHECK_FALSE(joints.is_valid(xc));
    CHECK(attached(joints, bodies, r2) == std::vector<sqp::JointHandle>{rx});
    CHECK(attached(joints, bodies, x1) == std::vector<sqp::JointHandle>{rx});
    REQUIRE(bodies.clear(joints));
    CHECK_FALSE(joints.is_valid(rx));
    CHECK(bodies.body_count() == 0);
    (void)r0;
}

TEST_CASE("Joint removal respects island preflight", "[joint][island]") {
    sqp::BodyPool bodies;
    sqp::IslandPool islands;
    sqp::JointPool joints;
    const auto a = bodies.add_rigid();
    const auto b = bodies.add_rigid();
    const auto island = islands.create();
    REQUIRE(islands.add_body(island, a, bodies));
    sqp::PinJointDef def;
    def.body_a = a;
    def.body_b = b;
    const auto joint = joints.create(bodies, def);
    REQUIRE(joint.valid());
    CHECK_FALSE(bodies.remove(a, joints));
    CHECK_FALSE(bodies.remove(a, islands));
    CHECK(joints.is_valid(joint));
    REQUIRE(bodies.remove(a, islands, joints));
    CHECK_FALSE(joints.is_valid(joint));
    CHECK_FALSE(bodies.is_valid(a));
    CHECK(bodies.is_valid(b));
    CHECK(bodies.clear(islands, joints));
}

TEST_CASE("JointPool refuses a different BodyPool with identical handle bits", "[joint][pool]") {
    sqp::BodyPool original;
    sqp::BodyPool other;
    sqp::JointPool joints;
    sqp::PinJointDef def;
    def.body_a = original.add_rigid();
    def.body_b = original.add_rigid();
    REQUIRE(other.add_rigid() == def.body_a);
    REQUIRE(other.add_rigid() == def.body_b);
    const auto joint = joints.create(original, def);
    REQUIRE(joint.valid());
    CHECK_FALSE(joints.create(other, def).valid());
    CHECK_FALSE(joints.destroy(other, joint));
    CHECK_FALSE(joints.first_attached(other, def.body_a).valid());
    CHECK(joints.is_valid(joint));
    REQUIRE(joints.destroy(original, joint));
    CHECK(joints.create(other, def).valid());
}
