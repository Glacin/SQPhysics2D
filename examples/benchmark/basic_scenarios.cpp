#include <limits>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "scenarios.h"
#include "scene_utils.h"
#include "world/world_stages.h"

namespace sqp::benchmark {
void BoxFreeMotionScenario::setup(Context& context, const Options& options) {
    if (options.body_count == 0)
        throw std::invalid_argument("no-collision body count must be positive");

    context.settings.gravity = {};
    context.settings.linear_velocity_tick_multiplier = 1;
    context.settings.angular_velocity_tick_multiplier = 1;
    bodies_.clear();
    bodies_.reserve(options.body_count);
    validated_ticks_ = 0;

    const auto square = scene_utils::rectangle(real(0.5), real(0.5));
    for (u32 i = 0; i < options.body_count; ++i) {
        RigidBodyDef box;
        box.position = {real(i) * real(3), 0};
        box.linear_velocity = {0, real(-5)};
        box.angular_velocity = 0;
        box.drag_multiplier = 1;
        bodies_.push_back(scene_utils::rigid_polygon(context, box, square));
    }
}

std::string BoxFreeMotionScenario::validate(const Context& context) const {
    ++validated_ticks_;
    const real expected_y = -real(5) * real(validated_ticks_) / real(config::TICKS_PER_SECOND);
    const real y_tolerance = real(0.001) + math::abs(expected_y) * real(0.0001);
    if (!context.state.contact_graph.contacts().empty())
        return "no-collision scene produced contacts";
    for (std::size_t i = 0; i < bodies_.size(); ++i) {
        BodyState state;
        if (!context.state.body_pool.get_state(bodies_[i], state) || !state.is_alive() ||
            !scene_utils::finite(state))
            return "no-collision box missing or non-finite";
        if (const auto id = context.state.body_pool.island_of(bodies_[i]);
            id && context.state.island_pool.is_sleeping(*id))
            return "no-collision box fell asleep";
        if (math::abs(state.position.x - real(i) * real(3)) > real(0.001) ||
            math::abs(state.position.y - expected_y) > y_tolerance || state.rotation != 0 ||
            state.linear_velocity.x != 0 ||
            math::abs(state.linear_velocity.y + real(5)) > real(0.0001) ||
            state.angular_velocity != 0)
            return "no-collision box deviated from uniform downward motion";
    }
    return {};
}

void BoxFloorSleepScenario::setup(Context& context, const Options& options) {
    if (!options.body_count)
        throw std::invalid_argument("floor sleep body count must be positive");
    if (options.initial_sleep && !options.allow_sleep)
        throw std::invalid_argument("initial sleep requires --sleep on");
    if (u64(options.body_count) + options.awake_bodies > std::numeric_limits<u32>::max())
        throw std::invalid_argument("floor sleep total body count exceeds u32 range");
    sleeping_count_ = options.body_count;
    awake_count_ = options.awake_bodies;
    context.settings.gravity = {0, -10};
    PhysicsMaterial material;
    material.base_restitution = sqp::real(0);
    scene_utils::install_material(context, material);
    const real width = real(options.body_count - 1) * real(1.5) + 1;
    AabbBoxBodyDef floor;
    floor.position = {0, real(-0.5)};
    floor_ = scene_utils::aabb_box(
        context, floor,
        {width / 2 + 1 + (awake_count_ ? 10 + real(awake_count_) * real(1.5) : 0), real(0.5)}, 0);
    bodies_.clear();
    bodies_.reserve(static_cast<std::size_t>(options.body_count) + options.awake_bodies);
    for (u32 i = 0; i < options.body_count; ++i) {
        RigidBodyDef box;
        box.position = {-real(options.body_count - 1) * real(0.75) + real(i) * real(1.5),
                        real(0.5)};
        bodies_.push_back(
            scene_utils::rigid_polygon(context, box, scene_utils::rectangle(real(0.5), real(0.5))));
        if (options.initial_sleep &&
            !park_island(context.state, *context.state.body_pool.island_of(bodies_.back())))
            throw std::runtime_error("initial sleep failed");
    }
    for (u32 i = 0; i < awake_count_; ++i) {
        RigidBodyDef box;
        box.position = {real(sleeping_count_ - 1) * real(0.75) + 10 + real(i) * real(1.5),
                        real(0.5)};
        box.configs.override_flags = BodyConfigs::OverrideFlags::NoSleep;
        bodies_.push_back(
            scene_utils::rigid_polygon(context, box, scene_utils::rectangle(real(0.5), real(0.5))));
    }
}

std::string BoxFloorSleepScenario::validate(const Context& context) const {
    BodyState floor;
    if (!context.state.body_pool.get_state(floor_, floor) || !floor.is_alive() ||
        !scene_utils::finite(floor))
        return "floor missing or non-finite";
    if (floor.position.x != 0 || floor.position.y != real(-0.5) || floor.rotation != 0 ||
        floor.linear_velocity.length() != 0 || floor.angular_velocity != 0)
        return "floor moved";
    const real tolerance = real(0.04);
    real previous_right = 0;
    for (std::size_t i = 0; i < bodies_.size(); ++i) {
        BodyState state;
        if (!context.state.body_pool.get_state(bodies_[i], state) || !state.is_alive() ||
            !scene_utils::finite(state))
            return "floor box missing or non-finite";
        const real initial_x = i < sleeping_count_
                                   ? -real(sleeping_count_ - 1) * real(0.75) + real(i) * real(1.5)
                                   : real(sleeping_count_ - 1) * real(0.75) + 10 +
                                         real(i - sleeping_count_) * real(1.5);
        if (math::abs(state.position.x - initial_x) > tolerance ||
            math::abs(state.position.y - real(0.5)) > tolerance)
            return "floor box drift exceeds 0.04 m";
        // Rotated square bounds, without allocating shape vectors during validation.
        const real extent = real(0.5) * (math::abs(math::cos(state.rotation)) +
                                         math::abs(math::sin(state.rotation)));
        if (math::abs(state.position.y - extent) > tolerance)
            return "floor box no longer rests on floor";
        if (i && state.position.x - extent <= previous_right)
            return "floor boxes touch or overlap";
        previous_right = state.position.x + extent;
        if (!context.settings.allow_sleep &&
            context.state.island_pool.is_sleeping(*context.state.body_pool.island_of(bodies_[i])))
            return "box asleep with sleep disabled";
    }
    return {};
}

SleepStats BoxFloorSleepScenario::sleep_stats(const Context& context) const {
    SleepStats stats;
    for (BodyHandle body : bodies_) {
        if (context.state.body_pool.is_valid(body)) {
            const auto island = context.state.body_pool.island_of(body);
            if (!island || !context.state.island_pool.is_sleeping(*island))
                ++stats.awake;
            else
                ++stats.sleeping;
        }
    }
    return stats;
}

std::string BoxFloorSleepScenario::validate_sleep(const Context& context,
                                                  tick_t elapsed_ticks) const {
    // Allow a settling window followed by a full quiet-time window.
    const tick_t threshold = context.settings.sleep_time_threshold_tu / config::TU_PER_TICK +
                             (context.settings.sleep_time_threshold_tu % config::TU_PER_TICK != 0);
    if (context.settings.allow_sleep && elapsed_ticks >= 2 * threshold + 1 &&
        sleep_stats(context).awake != awake_count_)
        return "floor boxes failed to sleep after settling and quiet-time windows";
    return {};
}
void PyramidScenario::setup_pyramid(Context& context, const Options& options, bool rigid_boxes) {
    const u32 layers = options.body_count;
    if (layers == 0)
        throw std::invalid_argument("pyramid bottom-row body count must be positive");
    const u64 body_count = u64(layers) * (u64(layers) + 1) / 2;
    if (body_count > std::numeric_limits<u32>::max())
        throw std::invalid_argument("pyramid body count exceeds u32 range");

    context.settings.gravity = {0, -10};
    PhysicsMaterial material;
    material.base_restitution = sqp::real(0);
    scene_utils::install_material(context, material);

    AabbBoxBodyDef platform_definition;
    platform_definition.position = {0, real(-0.5)};
    platform_ = scene_utils::aabb_box(context, platform_definition,
                                      {static_cast<real>(layers), real(0.5)}, 0);
    bodies_.clear();
    bodies_.reserve(static_cast<std::size_t>(body_count));
    layers_ = layers;
    top_body_ = {};
    initial_top_position_ = {};

    for (u32 row = 0; row < layers; ++row) {
        const u32 count = layers - row;
        const real start_x = -real(count - 1) * real(0.5);
        const real y = real(0.5) + static_cast<real>(row);
        for (u32 column = 0; column < count; ++column) {
            const Vec2 position{start_x + static_cast<real>(column), y};
            BodyHandle body;
            if (rigid_boxes) {
                RigidBodyDef definition;
                definition.position = position;
                body = scene_utils::rigid_polygon(context, definition,
                                                  scene_utils::rectangle(real(0.5), real(0.5)));
            } else {
                AabbBoxBodyDef definition;
                definition.position = position;
                body = scene_utils::aabb_box(context, definition, {real(0.5), real(0.5)});
            }
            bodies_.push_back(body);
            if (row + 1 == layers) {
                top_body_ = body;
                initial_top_position_ = position;
            }
        }
    }
}

PositionStats PyramidScenario::position_stats(const Context& context) const {
    BodyState state;
    if (!is_valid(top_body_) || !context.state.body_pool.get_state(top_body_, state))
        return {};
    PositionStats stats;
    stats.available = true;
    stats.position_delta = (state.position - initial_top_position_).length();
    stats.position_limit = static_cast<real>(layers_) * real(0.2);
    return stats;
}

std::vector<real> PyramidScenario::layer_horizontal_drift(const Context& context) const {
    std::vector<real> drift(layers_, real(0));
    std::size_t index = 0;
    for (u32 row = 0; row < layers_; ++row) {
        const u32 count = layers_ - row;
        const real start_x = -real(count - 1) * real(0.5);
        for (u32 column = 0; column < count; ++column) {
            BodyState state;
            if (!context.state.body_pool.get_state(bodies_[index++], state))
                continue;
            drift[row] = std::max(drift[row],
                                  std::abs(state.position.x - start_x - real(column)));
        }
    }
    return drift;
}

std::vector<real> PyramidScenario::layer_center_sag(const Context& context) const {
    std::vector<real> sag(layers_, real(0));
    std::size_t index = 0;
    for (u32 row = 0; row < layers_; ++row) {
        const u32 count = layers_ - row;
        for (u32 column = 0; column < count; ++column) {
            BodyState state;
            const bool valid = context.state.body_pool.get_state(bodies_[index++], state);
            if (!valid || (column != (count - 1) / 2 && column != count / 2))
                continue;
            sag[row] = std::max(sag[row], real(0.5) + real(row) - state.position.y);
        }
    }
    return sag;
}

std::string PyramidScenario::validate(const Context& context) const {
    BodyState platform;
    if (!context.state.body_pool.get_state(platform_, platform) || !platform.is_alive())
        return "pyramid platform missing";
    if (!scene_utils::finite(platform))
        return "non-finite pyramid platform state";
    for (const BodyHandle body : bodies_) {
        BodyState state;
        if (!context.state.body_pool.get_state(body, state) || !state.is_alive())
            return "pyramid body missing";
        if (!scene_utils::finite(state))
            return "non-finite pyramid body state";
    }
    const PositionStats stats = position_stats(context);
    if (!stats.available)
        return "pyramid top body missing";
    if (stats.position_delta > stats.position_limit) {
        std::ostringstream error;
        error << "top body moved too far: delta=" << stats.position_delta
              << " limit=" << stats.position_limit;
        return error.str();
    }
    return {};
}

} // namespace sqp::benchmark
