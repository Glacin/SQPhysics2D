#include <array>
#include <span>

#include "solver/contact_material.h"
#include "sqphysics2d/detail/contact_solver.h"
#include "sqphysics2d/detail/data_storage.h"
#include <catch_amalgamated.hpp>

namespace {
sqp::BodyHandle body(sqp::BodyPool& pool, sqp::ShapePool& shapes, sqp::Vec2 position,
                     sqp::Vec2 velocity,
                     sqp::BodyConfigs::RegisterType type = sqp::BodyConfigs::RegisterType::Dynamic,
                     sqp::u32 flags = 0) {
    const auto shape = shapes.add_circle({sqp::real(0.5), {}, 0, sqp::real(1)});
    sqp::RigidBodyDef def(position, 0, velocity, 0);
    def.configs.register_type = type;
    def.configs.override_flags = sqp::BodyConfigs::OverrideFlags(flags);
    const auto handle = pool.add_rigid(def);
    REQUIRE(pool.attach_shape(handle, shape, shapes));
    return handle;
}

sqp::BodyState state(const sqp::BodyPool& pool, sqp::BodyHandle handle) {
    sqp::BodyState result;
    REQUIRE(pool.get_state(handle, result));
    return result;
}

sqp::ToiContact contact(sqp::BodyHandle a, sqp::BodyHandle b, sqp::Vec2 point, sqp::Vec2 normal,
                        sqp::PhysicsMaterial material = {}) {
    sqp::ToiContact result;
    result.a = a;
    result.b = b;
    result.manifold.count = 1;
    result.manifold.normal = normal;
    result.manifold.points[0].point = point;
    result.manifold.points[0].depth = sqp::real(0.01);
    const auto combined = sqp::detail::combine_contact_materials(material, material);
    result.parameters = {material.base_restitution,     combined.dynamic_friction,
                         combined.static_friction,      combined.contact_damping_ratio,
                         combined.contact_frequency_hz, combined.max_penetration_correction_speed};
    return result;
}

std::vector<sqp::ContactSolveReport>
solve(sqp::BodyPool& pool, std::span<const sqp::ToiContact> contacts, sqp::u32 passes = 4) {
    std::vector<sqp::ContactSolveReport> reports;
    REQUIRE(sqp::solve_toi_contacts(pool, contacts, passes, sqp::real(0.1), reports));
    return reports;
}
} // namespace

TEST_CASE("TOI contact solver resolves rigid impacts and ignores separating contacts",
          "[contact][toi]") {
    sqp::BodyPool pool;
    sqp::ShapePool shapes;
    const auto wall = body(pool, shapes, {}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto moving = body(pool, shapes, {0, 1}, {0, -2});
    sqp::PhysicsMaterial material;
    material.base_restitution = sqp::real(1);
    material.restitution_flags = sqp::RESTITUTION_DISABLE_IMPACT_SCALING;
    const auto hit = contact(wall, moving, {0, sqp::real(0.5)}, {0, 1}, material);
    const auto before_pose = state(pool, moving);
    const auto reports = solve(pool, std::span<const sqp::ToiContact>(&hit, 1));
    CHECK(state(pool, moving).position == before_pose.position);
    CHECK(state(pool, moving).rotation == before_pose.rotation);
    CHECK(state(pool, moving).linear_velocity.y == Catch::Approx(2));
    REQUIRE(reports.size() == 1);
    CHECK(reports[0].restitution_applied);
    CHECK(reports[0].normal_impulse > 0);
    const auto separated = solve(pool, std::span<const sqp::ToiContact>(&hit, 1));
    CHECK(state(pool, moving).linear_velocity.y == Catch::Approx(2));
    CHECK(separated[0].normal_impulse == 0);
}

TEST_CASE("TOI contact solver transfers momentum between dynamic bodies", "[contact][toi]") {
    sqp::BodyPool pool;
    sqp::ShapePool shapes;
    const auto a = body(pool, shapes, {-sqp::real(0.5), 0}, {2, 0});
    const auto b = body(pool, shapes, {sqp::real(0.5), 0}, {});
    sqp::PhysicsMaterial material;
    material.base_restitution = sqp::real(1);
    material.restitution_flags = sqp::RESTITUTION_DISABLE_IMPACT_SCALING;
    const auto hit = contact(a, b, {}, {1, 0}, material);
    solve(pool, std::span<const sqp::ToiContact>(&hit, 1));
    CHECK(state(pool, a).linear_velocity.x == Catch::Approx(0).margin(1e-5));
    CHECK(state(pool, b).linear_velocity.x == Catch::Approx(2).margin(1e-5));
}

TEST_CASE("TOI contact solver applies friction and angular response while respecting locks",
          "[contact][toi]") {
    bool responses_valid = true;
    for (bool lock_rotation : {false, true}) {
        sqp::BodyPool pool;
        sqp::ShapePool shapes;
        const auto wall =
            body(pool, shapes, {}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
        const auto moving = body(
            pool, shapes, {0, 1}, {2, -2}, sqp::BodyConfigs::RegisterType::Dynamic,
            lock_rotation ? static_cast<sqp::u32>(sqp::BodyConfigs::OverrideFlags::LockRot) : 0);
        sqp::PhysicsMaterial material;
        material.static_friction = sqp::real(1);
        material.dynamic_friction = sqp::real(1);
        const auto hit = contact(wall, moving, {sqp::real(0.4), sqp::real(0.5)}, {0, 1}, material);
        const auto reports = solve(pool, std::span<const sqp::ToiContact>(&hit, 1));
        responses_valid &= state(pool, moving).linear_velocity.x < 2 && !reports.empty() &&
                           reports[0].tangent_impulse > 0 &&
                           (lock_rotation ? state(pool, moving).angular_velocity == 0
                                          : state(pool, moving).angular_velocity != 0);
    }
    CHECK(responses_valid);
}

TEST_CASE("TOI contact solver does not impulse a locked normal axis", "[contact][toi]") {
    sqp::BodyPool pool;
    sqp::ShapePool shapes;
    const auto wall = body(pool, shapes, {}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto moving = body(pool, shapes, {0, 1}, {0, -2}, sqp::BodyConfigs::RegisterType::Dynamic,
                             static_cast<sqp::u32>(sqp::BodyConfigs::OverrideFlags::LockLinearY));
    const auto hit = contact(wall, moving, {0, sqp::real(0.5)}, {0, 1});
    const auto reports = solve(pool, std::span<const sqp::ToiContact>(&hit, 1));
    CHECK(state(pool, moving).linear_velocity.y == 0);
    CHECK(reports[0].normal_impulse == 0);
}

TEST_CASE("TOI contact solver handles two points and a simultaneous contact chain",
          "[contact][toi]") {
    sqp::BodyPool pool;
    sqp::ShapePool shapes;
    const auto wall = body(pool, shapes, {}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto a = body(pool, shapes, {0, 1}, {0, -3});
    const auto b = body(pool, shapes, {0, 2}, {0, -3});
    auto first = contact(wall, a, {-sqp::real(0.25), sqp::real(0.5)}, {0, 1});
    first.manifold.count = 2;
    first.manifold.points[1] = first.manifold.points[0];
    first.manifold.points[1].point.x = sqp::real(0.25);
    const auto second = contact(a, b, {0, sqp::real(1.5)}, {0, 1});
    const std::array contacts{first, second};
    const auto reports = solve(pool, contacts, 8);
    REQUIRE(reports.size() == 2);
    CHECK(reports[0].total_impulse >= reports[0].normal_impulse);
    CHECK(reports[0].normal_impulse > 0);
    CHECK(state(pool, a).linear_velocity.y >= -sqp::real(0.01));
    CHECK(state(pool, b).linear_velocity.y >= -sqp::real(0.01));
}

TEST_CASE("TOI contact solver rejects stale and reused body handles",
          "[contact][toi][validation]") {
    sqp::BodyPool pool;
    sqp::ShapePool shapes;
    const auto wall = body(pool, shapes, {}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto removed = body(pool, shapes, {0, 1}, {0, -1});
    const auto hit = contact(wall, removed, {0, sqp::real(0.5)}, {0, 1});
    REQUIRE(pool.remove(removed));
    const auto replacement = body(pool, shapes, {2, 1}, {});
    REQUIRE(replacement.index() == removed.index());
    std::vector<sqp::ContactSolveReport> reports;
    CHECK_FALSE(sqp::solve_toi_contacts(pool, std::span<const sqp::ToiContact>(&hit, 1), 1,
                                        sqp::real(0.1), reports));
    CHECK(reports.empty());
    CHECK(state(pool, replacement).linear_velocity == sqp::Vec2{});
}

TEST_CASE("TOI contact solver skips disabled zero-point contacts", "[contact][toi]") {
    sqp::BodyPool pool;
    sqp::ShapePool shapes;
    const auto wall = body(pool, shapes, {}, {}, sqp::BodyConfigs::RegisterType::ImmutableStatic);
    const auto moving = body(pool, shapes, {0, 1}, {0, -2});
    auto hit = contact(wall, moving, {0, sqp::real(0.5)}, {0, 1});
    hit.manifold.count = 0;
    const auto reports = solve(pool, std::span<const sqp::ToiContact>(&hit, 1));
    CHECK(reports.empty());
    CHECK(state(pool, moving).linear_velocity.y == -2);
}
