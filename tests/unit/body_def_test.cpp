#include "sqphysics2d/body.h"
#include <catch_amalgamated.hpp>

TEST_CASE("body definitions preserve defaults and constructor inputs", "[body][definition]") {
    const sqp::RigidBodyDef default_rigid;
    CHECK(default_rigid.rotation == 0);
    CHECK(default_rigid.angular_velocity == 0);
    CHECK(default_rigid.drag_multiplier == 1);
    CHECK(default_rigid.configs.register_type == sqp::BodyConfigs::RegisterType::Dynamic);
    CHECK(default_rigid.call_back_ptr == nullptr);

    const sqp::RigidBodyDef rigid({1, 2}, 3, {4, 5}, 6);
    CHECK(rigid.position == sqp::Vec2(1, 2));
    CHECK(rigid.rotation == 3);
    CHECK(rigid.linear_velocity == sqp::Vec2(4, 5));
    CHECK(rigid.angular_velocity == 6);

    const sqp::LiteCircleBodyDef circle({1, 2}, {3, 4}, 7);
    CHECK(circle.position == sqp::Vec2(1, 2));
    CHECK(circle.linear_velocity == sqp::Vec2(3, 4));
    CHECK(circle.drag_multiplier == 7);

    const sqp::AabbBoxBodyDef box({1, 2}, {3, 4}, 8);
    CHECK(box.position == sqp::Vec2(1, 2));
    CHECK(box.linear_velocity == sqp::Vec2(3, 4));
    CHECK(box.drag_multiplier == 8);
}

TEST_CASE("body configs provide typed defaults and combinable flags", "[body][configs]") {
    sqp::BodyConfigs configs;
    CHECK(configs.kind == sqp::BodyConfigs::Kind::None);
    CHECK(configs.register_type == sqp::BodyConfigs::RegisterType::Dynamic);
    CHECK(configs.ccd_override == sqp::BodyConfigs::CcdOverride::Auto);
    CHECK(configs.runtime_flags == sqp::BodyConfigs::RuntimeFlags(0));
    CHECK(configs.event_flags == sqp::BodyConfigs::EventFlags(0));
    CHECK(configs.override_flags == sqp::BodyConfigs::OverrideFlags(0));

    configs.event_flags =
        sqp::BodyConfigs::EventFlags::ContactEvents | sqp::BodyConfigs::EventFlags::AfterSolve;
    CHECK((configs.event_flags & sqp::BodyConfigs::EventFlags::ContactEvents) !=
          sqp::BodyConfigs::EventFlags(0));
    CHECK((configs.event_flags & sqp::BodyConfigs::EventFlags::SleepEvents) ==
          sqp::BodyConfigs::EventFlags(0));
    configs.override_flags =
        sqp::BodyConfigs::OverrideFlags::LockLinearX | sqp::BodyConfigs::OverrideFlags::NoGravityY;
    CHECK((configs.override_flags & sqp::BodyConfigs::OverrideFlags::LockLinearX) !=
          sqp::BodyConfigs::OverrideFlags(0));
}

TEST_CASE("collision filters require mutual category acceptance", "[body][filter]") {
    constexpr sqp::u64 low = sqp::u64(1) << 1;
    constexpr sqp::u64 high = sqp::u64(1) << 7;
    constexpr sqp::CollisionFilter a{low, high};
    constexpr sqp::CollisionFilter b{high, low};
    constexpr sqp::CollisionFilter rejects_a{high, high};
    STATIC_REQUIRE(sqp::should_collide(a, b));
    STATIC_REQUIRE_FALSE(sqp::should_collide(a, rejects_a));
    STATIC_REQUIRE_FALSE(sqp::should_collide(rejects_a, a));
}

TEST_CASE("body override utilities apply typed policies", "[body][flags][ccd]") {
    CHECK(sqp::gravity_with_body_overrides({2, 3}, sqp::BodyConfigs::OverrideFlags(0)) ==
          sqp::Vec2(2, 3));
    CHECK(sqp::gravity_with_body_overrides({2, 3}, sqp::BodyConfigs::OverrideFlags::NoGravityX) ==
          sqp::Vec2(0, 3));
    CHECK(sqp::gravity_with_body_overrides({2, 3}, sqp::BodyConfigs::OverrideFlags::NoGravityY) ==
          sqp::Vec2(2, 0));
    CHECK(sqp::gravity_with_body_overrides({2, 3}, sqp::BodyConfigs::OverrideFlags::NoGravity) ==
          sqp::Vec2());

    bool overrides_preserved = true;
    for (const auto type :
         {sqp::BodyConfigs::CcdOverride::Auto, sqp::BodyConfigs::CcdOverride::Shapecast,
          sqp::BodyConfigs::CcdOverride::Disabled}) {
        sqp::BodyConfigs value;
        value.ccd_override = type;
        overrides_preserved &= value.ccd_override == type;
    }
    CHECK(overrides_preserved);
    CHECK(static_cast<sqp::u8>(sqp::BodyConfigs::CcdOverride::Auto) == 0);
    CHECK(static_cast<sqp::u8>(sqp::BodyConfigs::CcdOverride::Shapecast) == 1);
    CHECK(static_cast<sqp::u8>(sqp::BodyConfigs::CcdOverride::Disabled) == 2);
}

TEST_CASE("lite body definitions only describe body state and policies", "[body][definition]") {
    const sqp::LiteCircleBodyDef circle;
    const sqp::AabbBoxBodyDef box;
    CHECK(circle.position == sqp::Vec2());
    CHECK(box.position == sqp::Vec2());
    CHECK(circle.linear_velocity == sqp::Vec2());
    CHECK(box.linear_velocity == sqp::Vec2());
    CHECK(circle.drag_multiplier == 1);
    CHECK(box.drag_multiplier == 1);
    CHECK(circle.configs.override_flags == sqp::BodyConfigs::OverrideFlags(0));
    CHECK(box.configs.override_flags == sqp::BodyConfigs::OverrideFlags(0));
}
