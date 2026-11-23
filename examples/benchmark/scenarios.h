#pragma once

#include <memory>
#include <string>

#include "scenario_runner.h"

namespace sqp::benchmark {

class BoxFloorSleepScenario final : public Scenario {
public:
    const char* name() const override {
        return "box_floor_sleep";
    }
    void setup(Context&, const Options&) override;
    std::string validate(const Context&) const override;
    SleepStats sleep_stats(const Context&) const;
    std::string validate_sleep(const Context&, tick_t) const override;

private:
    std::vector<BodyHandle> bodies_;
    BodyHandle floor_;
    u32 sleeping_count_ = 0, awake_count_ = 0;
};

class PyramidScenario : public Scenario {
public:
    std::string validate(const Context&) const override;
    PositionStats position_stats(const Context&) const;
    std::vector<real> layer_horizontal_drift(const Context&) const override;
    std::vector<real> layer_center_sag(const Context&) const override;

protected:
    void setup_pyramid(Context&, const Options&, bool rigid_boxes);

private:
    std::vector<BodyHandle> bodies_;
    BodyHandle platform_;
    BodyHandle top_body_;
    Vec2 initial_top_position_;
    u32 layers_ = 0;
};

class AabbBoxPyramidScenario final : public PyramidScenario {
public:
    const char* name() const override {
        return "aabb_box_pyramid";
    }
    void setup(Context& context, const Options& options) override {
        setup_pyramid(context, options, false);
    }
};

class RigidBoxPyramidScenario final : public PyramidScenario {
public:
    const char* name() const override {
        return "rigid_box_pyramid";
    }
    void setup(Context& context, const Options& options) override {
        setup_pyramid(context, options, true);
    }
};

class WallBounceScenario : public Scenario {
public:
    std::string validate(const Context&) const override;

protected:
    struct Record {
        BodyHandle body{};
        real radius = 0;
        mutable Vec2 previous_position;
    };
    void track_body(const Context&, BodyHandle, real radius);
    std::vector<Record> bodies_;
    real validation_extent_ = 0;
};

class CircleBoxBounceScenario final : public WallBounceScenario {
public:
    const char* name() const override {
        return "circle_box_bounce";
    }
    void setup(Context&, const Options&) override;
};

class ConvexPolygonBounceScenario final : public WallBounceScenario {
public:
    const char* name() const override {
        return "convex_polygon_bounce";
    }
    void setup(Context&, const Options&) override;
};

class MixedShapeBounceScenario final : public WallBounceScenario {
public:
    const char* name() const override {
        return "mixed_shape_bounce";
    }
    void setup(Context&, const Options&) override;
};

class BoxFreeMotionScenario final : public Scenario {
public:
    const char* name() const override {
        return "box_free_motion";
    }
    void setup(Context&, const Options&) override;
    std::string validate(const Context&) const override;
    bool default_allow_sleep() const override {
        return false;
    }

private:
    std::vector<BodyHandle> bodies_;
    mutable tick_t validated_ticks_ = 0;
};

std::unique_ptr<Scenario> make_scenario(const std::string& name);

} // namespace sqp::benchmark
