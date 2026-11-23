#include "sqphysics2d/world.h"

#include <memory>
#include <type_traits>
#include <utility>

#include <catch_amalgamated.hpp>

namespace {
static_assert(sizeof(sqp::BodyPool) > 0);
static_assert(sizeof(sqp::ShapePool) > 0);

template <class T> void check_command_alternative() {
    sqp::PhysicsCommand command{std::in_place_type<T>};
    CHECK(std::holds_alternative<T>(command));
}

template <class T> void check_body_definition() {
    auto definition = std::make_unique<T>();
    const auto* original = definition.get();
    sqp::CreateBody create{std::move(definition)};
    sqp::PhysicsCommand command{std::move(create)};
    const auto& stored = std::get<sqp::CreateBody>(command);
    const auto& owned = std::get<std::unique_ptr<T>>(stored.definition);
    REQUIRE(owned);
    CHECK(owned.get() == original);
}
} // namespace

TEST_CASE("PhysicsCommand accepts every built-in command", "[world][command]") {
    check_command_alternative<sqp::SetVelocity>();
    check_command_alternative<sqp::SetPose>();
    check_command_alternative<sqp::SetBodyConfig>();
    check_command_alternative<sqp::RemoveBody>();
    check_command_alternative<sqp::ApplyImpulse>();
    check_command_alternative<sqp::SetConstantForce>();
    check_command_alternative<sqp::ClearForces>();
    check_command_alternative<sqp::SetCollisionFilter>();
    check_command_alternative<sqp::AttachShape>();
    check_command_alternative<sqp::DetachShape>();
    check_command_alternative<sqp::CreateBody>();
    check_command_alternative<sqp::CreateShape>();
    check_command_alternative<sqp::RemoveShape>();
    check_command_alternative<sqp::SetSettings>();
    check_command_alternative<sqp::SetCallback>();
}

TEST_CASE("SetBodyConfig carries each supported configuration change", "[world][command]") {
    using Config = sqp::SetBodyConfig;
    const sqp::BodyHandle body{7, 2, sqp::BodyConfigs::Kind::Rigid};
    Config command{body, sqp::BodyConfigs::CcdOverride::Auto};
    CHECK(command.body == body);
    command.value = sqp::BodyConfigs::CcdOverride::Shapecast;
    CHECK(std::get<sqp::BodyConfigs::CcdOverride>(command.value) ==
          sqp::BodyConfigs::CcdOverride::Shapecast);
    command.value = sqp::BodyConfigs::EventFlags::ContactEvents;
    CHECK(std::get<sqp::BodyConfigs::EventFlags>(command.value) ==
          sqp::BodyConfigs::EventFlags::ContactEvents);
    command.value = sqp::BodyConfigs::OverrideFlags::NoGravity;
    CHECK(std::get<sqp::BodyConfigs::OverrideFlags>(command.value) ==
          sqp::BodyConfigs::OverrideFlags::NoGravity);
    sqp::PhysicsCommand stored{command};
    CHECK(std::holds_alternative<Config>(stored));
}

TEST_CASE("CreateBody owns exactly one body definition", "[world][command]") {
    static_assert(!std::is_copy_constructible_v<sqp::CreateBody>);
    static_assert(std::is_move_constructible_v<sqp::CreateBody>);
    check_body_definition<sqp::RigidBodyDef>();
    check_body_definition<sqp::LiteCircleBodyDef>();
    check_body_definition<sqp::AabbBoxBodyDef>();
}
