#include "scenarios.h"

namespace sqp::benchmark {

std::unique_ptr<Scenario> make_scenario(const std::string& name) {
    if (name == "box_floor_sleep" || name == "rigid_box_floor_sleep")
        return std::make_unique<BoxFloorSleepScenario>();
    if (name == "circle_box_bounce" || name == "lite_wall_bounce")
        return std::make_unique<CircleBoxBounceScenario>();
    if (name == "convex_polygon_bounce" || name == "polygon_wall_bounce")
        return std::make_unique<ConvexPolygonBounceScenario>();
    if (name == "mixed_shape_bounce" || name == "mixed_wall_bounce")
        return std::make_unique<MixedShapeBounceScenario>();
    if (name == "aabb_box_pyramid" || name == "aabb_pyramid")
        return std::make_unique<AabbBoxPyramidScenario>();
    if (name == "rigid_box_pyramid")
        return std::make_unique<RigidBoxPyramidScenario>();
    if (name == "box_free_motion")
        return std::make_unique<BoxFreeMotionScenario>();
    return {};
}

} // namespace sqp::benchmark
