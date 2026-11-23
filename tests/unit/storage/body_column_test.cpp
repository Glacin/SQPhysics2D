#include <array>
#include <type_traits>
#include <utility>

#include "sqphysics2d/detail/data_storage.h"
#include <catch_amalgamated.hpp>

using namespace sqp;

static_assert(std::is_same_v<decltype(std::declval<BodyPool&>().scalar_column(
                                 BodyConfigs::Kind::Rigid, BodyScalarColumn::PositionX)),
                             std::span<real>>);
static_assert(std::is_same_v<decltype(std::declval<const BodyPool&>().scalar_column(
                                 BodyConfigs::Kind::Rigid, BodyScalarColumn::PositionX)),
                             std::span<const real>>);
static_assert(std::is_same_v<
              decltype(std::declval<const BodyPool&>().configs_column(BodyConfigs::Kind::Rigid)),
              std::span<const BodyConfigs>>);
static_assert(std::is_same_v<
              decltype(std::declval<BodyPool&>().motion_version_column(BodyConfigs::Kind::Rigid)),
              std::span<u64>>);
static_assert(std::is_same_v<decltype(std::declval<const BodyPool&>().motion_version_column(
                                 BodyConfigs::Kind::Rigid)),
                             std::span<const u64>>);

TEST_CASE("BodyPool exposes every stored scalar column by kind", "[body][columns]") {
    BodyPool bodies;
    ShapePool shapes;
    CircleShapeDef circle;
    circle.radius = real(0.5);
    PolygonShapeDef box;
    box.vertices = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
    const auto rigid = bodies.add_rigid({});
    const auto lite = bodies.add_lite_circle({}, shapes.add_circle(circle), shapes);
    const auto aabb = bodies.add_aabb_box({}, shapes.add_polygon(box), shapes);
    REQUIRE(bodies.is_valid(rigid));
    REQUIRE(bodies.is_valid(lite));
    REQUIRE(bodies.is_valid(aabb));

    constexpr std::array common{BodyScalarColumn::PositionX,     BodyScalarColumn::PositionY,
                                BodyScalarColumn::VelocityX,     BodyScalarColumn::VelocityY,
                                BodyScalarColumn::Mass,          BodyScalarColumn::InverseMass,
                                BodyScalarColumn::DragMultiplier};
    constexpr std::array rigid_only{
        BodyScalarColumn::Rotation,     BodyScalarColumn::AngularVelocity,
        BodyScalarColumn::Inertia,      BodyScalarColumn::InverseInertia,
        BodyScalarColumn::LocalCenterX, BodyScalarColumn::LocalCenterY,
        BodyScalarColumn::BoundRadius};
    constexpr std::array lite_only{BodyScalarColumn::Density, BodyScalarColumn::Radius,
                                   BodyScalarColumn::HalfExtentX, BodyScalarColumn::HalfExtentY};
    for (auto kind :
         {BodyConfigs::Kind::Rigid, BodyConfigs::Kind::LiteCircle, BodyConfigs::Kind::AabbBox}) {
        for (auto field : common) {
            CHECK(bodies.scalar_column(kind, field).size() == 1);
            CHECK(std::as_const(bodies).scalar_column(kind, field).size() == 1);
        }
        for (auto field : rigid_only)
            CHECK(bodies.scalar_column(kind, field).size() ==
                  (kind == BodyConfigs::Kind::Rigid ? 1 : 0));
        for (auto field : lite_only)
            CHECK(bodies.scalar_column(kind, field).size() ==
                  (kind == BodyConfigs::Kind::Rigid ? 0 : 1));
        CHECK(bodies.configs_column(kind).size() == 1);
        CHECK(bodies.motion_version_column(kind).size() == 1);
        CHECK(std::as_const(bodies).configs_column(kind).size() == 1);
        CHECK(std::as_const(bodies).motion_version_column(kind)[0] == 0);
    }
    CHECK(bodies.scalar_column(BodyConfigs::Kind::None, BodyScalarColumn::PositionX).empty());
    CHECK(
        bodies.scalar_column(BodyConfigs::Kind::Rigid, static_cast<BodyScalarColumn>(255)).empty());
    CHECK(std::as_const(bodies)
              .scalar_column(BodyConfigs::Kind::None, BodyScalarColumn::Mass)
              .empty());
    CHECK(bodies.configs_column(BodyConfigs::Kind::None).empty());
    CHECK(bodies.motion_version_column(BodyConfigs::Kind::None).empty());
}

TEST_CASE("BodyPool columns are writable", "[body][columns]") {
    BodyPool bodies;
    const auto body = bodies.add_rigid({});
    const auto index = bodies.dense_index(body);
    REQUIRE(index == 0);
    auto x = bodies.scalar_column(BodyConfigs::Kind::Rigid, BodyScalarColumn::PositionX);
    auto velocity = bodies.scalar_column(BodyConfigs::Kind::Rigid, BodyScalarColumn::VelocityX);
    auto configs = bodies.configs_column(BodyConfigs::Kind::Rigid);
    auto versions = bodies.motion_version_column(BodyConfigs::Kind::Rigid);
    x[*index] = real(7);
    velocity[*index] = real(3);
    configs[*index].override_flags |= BodyConfigs::OverrideFlags::LockLinearX;
    versions[*index] = u64{1} << 40;
    CHECK(std::as_const(bodies).scalar_column(BodyConfigs::Kind::Rigid,
                                              BodyScalarColumn::PositionX)[0] == 7);
    CHECK(std::as_const(bodies).scalar_column(BodyConfigs::Kind::Rigid,
                                              BodyScalarColumn::VelocityX)[0] == 3);
    CHECK((std::as_const(bodies).configs_column(BodyConfigs::Kind::Rigid)[0].override_flags &
           BodyConfigs::OverrideFlags::LockLinearX) != BodyConfigs::OverrideFlags(0));
    CHECK(std::as_const(bodies).motion_version_column(BodyConfigs::Kind::Rigid)[0] ==
          (u64{1} << 40));
}

TEST_CASE("BodyPool dense indices track removal and reject stale handles", "[body][columns]") {
    BodyPool bodies;
    RigidBodyDef first_def, second_def, third_def;
    first_def.position.x = 1;
    second_def.position.x = 2;
    third_def.position.x = 3;
    third_def.configs.ccd_override = BodyConfigs::CcdOverride::Shapecast;
    third_def.configs.event_flags = BodyConfigs::EventFlags::PreSolve;
    const auto first = bodies.add_rigid(first_def);
    const auto second = bodies.add_rigid(second_def);
    const auto third = bodies.add_rigid(third_def);
    REQUIRE(bodies.dense_index(first) == 0);
    REQUIRE(bodies.dense_index(second) == 1);
    REQUIRE(bodies.dense_index(third) == 2);
    bodies.motion_version_column(BodyConfigs::Kind::Rigid)[2] = 17;
    REQUIRE(bodies.remove(second));
    CHECK_FALSE(bodies.dense_index(second).has_value());
    REQUIRE(bodies.dense_index(third) == 1);
    CHECK(bodies.scalar_column(BodyConfigs::Kind::Rigid, BodyScalarColumn::PositionX)[1] == 3);
    CHECK(bodies.motion_version_column(BodyConfigs::Kind::Rigid)[1] == 17);
    CHECK(bodies.configs_column(BodyConfigs::Kind::Rigid)[1].ccd_override ==
          BodyConfigs::CcdOverride::Shapecast);
    CHECK(bodies.configs_column(BodyConfigs::Kind::Rigid)[1].event_flags ==
          BodyConfigs::EventFlags::PreSolve);
    const auto replacement = bodies.add_rigid({});
    REQUIRE(replacement.index() == second.index());
    CHECK_FALSE(bodies.dense_index(second).has_value());
    REQUIRE(bodies.dense_index(replacement) == 2);
    CHECK(bodies.motion_version_column(BodyConfigs::Kind::Rigid)[2] == 0);
    CHECK(bodies.configs_column(BodyConfigs::Kind::Rigid)[2].ccd_override ==
          BodyConfigs::CcdOverride::Auto);
    CHECK_FALSE(bodies.dense_index({}).has_value());
    bodies.clear();
    CHECK_FALSE(bodies.dense_index(first).has_value());
    CHECK_FALSE(bodies.dense_index(third).has_value());
    CHECK_FALSE(bodies.dense_index(replacement).has_value());
    CHECK(bodies.scalar_column(BodyConfigs::Kind::Rigid, BodyScalarColumn::PositionX).empty());
    CHECK(bodies.motion_version_column(BodyConfigs::Kind::Rigid).empty());
}

TEST_CASE("BodyPool columns and indices can be reacquired after move", "[body][columns]") {
    BodyPool source;
    RigidBodyDef def;
    def.position.x = real(9);
    const auto body = source.add_rigid(def);
    source.motion_version_column(BodyConfigs::Kind::Rigid)[0] = 23;
    BodyPool target(std::move(source));
    CHECK_FALSE(source.dense_index(body).has_value());
    CHECK(source.scalar_column(BodyConfigs::Kind::Rigid, BodyScalarColumn::PositionX).empty());
    CHECK(source.configs_column(BodyConfigs::Kind::Rigid).empty());
    CHECK(source.motion_version_column(BodyConfigs::Kind::Rigid).empty());
    REQUIRE(target.dense_index(body) == 0);
    CHECK(target.scalar_column(BodyConfigs::Kind::Rigid, BodyScalarColumn::PositionX)[0] == 9);
    CHECK(target.motion_version_column(BodyConfigs::Kind::Rigid)[0] == 23);
}
