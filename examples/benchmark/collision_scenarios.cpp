#include <algorithm>
#include <array>
#include <random>
#include <sstream>
#include <utility>
#include <vector>

#include "scenarios.h"
#include "scene_utils.h"

namespace sqp::benchmark {
namespace {

PhysicsMaterial bounce_material() {
    PhysicsMaterial material;
    material.base_restitution = sqp::real(1);
    material.restitution_flags = sqp::RESTITUTION_DISABLE_IMPACT_SCALING;
    material.dynamic_friction = sqp::real(0);
    material.static_friction = sqp::real(0);
    return material;
}

void lite_walls(Context& context, real half, real thickness) {
    for (const Vec2 position : std::array<Vec2, 4>{{{0, -half - thickness / 2},
                                                    {0, half + thickness / 2},
                                                    {-half - thickness / 2, 0},
                                                    {half + thickness / 2, 0}}}) {
        AabbBoxBodyDef definition;
        definition.position = position;
        const Vec2 extents =
            position.x == 0 ? Vec2(half + thickness, thickness / 2) : Vec2(thickness / 2, half);
        scene_utils::aabb_box(context, definition, extents, 0);
    }
}

void rigid_walls(Context& context, real half, real thickness) {
    const std::array<std::pair<Vec2, Vec2>, 4> walls = {
        {{{0, -half - thickness / 2}, {half + thickness, thickness / 2}},
         {{0, half + thickness / 2}, {half + thickness, thickness / 2}},
         {{-half - thickness / 2, 0}, {thickness / 2, half}},
         {{half + thickness / 2, 0}, {thickness / 2, half}}}};
    for (const auto& [position, extents] : walls) {
        RigidBodyDef definition;
        definition.position = position;
        definition.configs.register_type = BodyConfigs::RegisterType::ImmutableStatic;
        scene_utils::rigid_polygon(context, definition,
                                   scene_utils::rectangle(extents.x, extents.y));
    }
}

std::vector<Vec2> polygon(std::mt19937& rng, real min_area = 1, real max_area = 100) {
    std::uniform_int_distribution<int> side_count(3, 7);
    std::uniform_real_distribution<float> area(static_cast<float>(min_area),
                                               static_cast<float>(max_area));
    std::uniform_real_distribution<float> aspect(0.65f, 1.55f);
    std::uniform_real_distribution<float> offset(0, static_cast<float>(2 * math::PI));
    const int count = side_count(rng);
    const real target = area(rng), ratio = aspect(rng), angle_offset = offset(rng);
    std::vector<Vec2> result;
    result.reserve(static_cast<std::size_t>(count));
    for (int index = 0; index < count; ++index) {
        const real angle = 2 * math::PI * real(index) / real(count) + angle_offset;
        result.push_back({math::cos(angle) * ratio, math::sin(angle) / ratio});
    }
    real signed_area = 0;
    for (int index = 0; index < count; ++index)
        signed_area += math::cross(result[index], result[(index + 1) % count]);
    if (signed_area < 0) {
        std::reverse(result.begin(), result.end());
        signed_area = -signed_area;
    }
    const real scale = math::sqrt(target / (signed_area * real(0.5)));
    for (Vec2& point : result)
        point *= scale;
    return result;
}

real radius_of(const std::vector<Vec2>& vertices) {
    real radius = 0;
    for (const Vec2 point : vertices)
        radius = std::max(radius, point.length());
    return radius;
}

struct PlacedBody {
    Vec2 position;
    real radius;
};

Vec2 non_overlapping_position(std::mt19937& rng, real extent, real radius,
                              const std::vector<PlacedBody>& placed) {
    std::uniform_real_distribution<float> coordinate(static_cast<float>(-extent + radius),
                                                     static_cast<float>(extent - radius));
    Vec2 candidate;
    for (u32 attempt = 0; attempt < 64; ++attempt) {
        candidate = {coordinate(rng), coordinate(rng)};
        bool overlaps = false;
        for (const PlacedBody& other : placed) {
            const real separation = radius + other.radius;
            if ((candidate - other.position).length_squared() < separation * separation) {
                overlaps = true;
                break;
            }
        }
        if (!overlaps)
            return candidate;
    }
    return candidate;
}

Vec2 velocity(std::mt19937& rng, real low, real high) {
    std::uniform_real_distribution<float> speed(static_cast<float>(low), static_cast<float>(high));
    std::uniform_real_distribution<float> angle(0, static_cast<float>(2 * math::PI));
    const real magnitude = speed(rng), radians = angle(rng);
    return {math::cos(radians) * magnitude, math::sin(radians) * magnitude};
}

} // namespace

void WallBounceScenario::track_body(const Context& context, BodyHandle body, real radius) {
    BodyState state;
    if (!context.state.body_pool.get_state(body, state))
        throw std::runtime_error("bounce body missing during setup");
    bodies_.push_back({body, radius, state.position});
}

std::string WallBounceScenario::validate(const Context& context) const {
    for (const Record& record : bodies_) {
        BodyState state;
        if (!context.state.body_pool.get_state(record.body, state))
            return "body missing";
        if (!scene_utils::finite(state))
            return "non-finite body state";
        const auto island = context.state.body_pool.island_of(record.body);
        const bool frozen_this_tick =
            island && context.state.island_pool.last_frozen_until(*island) ==
                          context.state.schedule.time;
        if (state.position == record.previous_position && !frozen_this_tick) {
            std::ostringstream error;
            error << "body position did not update: handle=("
                  << static_cast<u32>(record.body.kind()) << ',' << record.body.index() << ','
                  << record.body.generation << ") position=(" << state.position.x << ','
                  << state.position.y << ')';
            return error.str();
        }
        record.previous_position = state.position;
        const Aabb bounds = scene_utils::body_bounds(context, record.body);
        if (!bounds.is_valid())
            return "body bounds missing";
        // The current solver can allow transient wall penetration; reject a body
        // once its whole bounds leave the enclosing walls.
        if (bounds.max.x < -validation_extent_ || bounds.min.x > validation_extent_ ||
            bounds.max.y < -validation_extent_ || bounds.min.y > validation_extent_) {
            std::ostringstream error;
            error << "body escaped: handle=(" << static_cast<u32>(record.body.kind()) << ','
                  << record.body.index() << ',' << record.body.generation << ") bounds=("
                  << bounds.min.x << ',' << bounds.min.y << ")-(" << bounds.max.x << ','
                  << bounds.max.y << ") position=(" << state.position.x << ',' << state.position.y
                  << ") velocity=(" << state.linear_velocity.x << ',' << state.linear_velocity.y
                  << ')';
            return error.str();
        }
    }
    return {};
}

void CircleBoxBounceScenario::setup(Context& context, const Options& options) {
    context.settings.gravity = {};
    scene_utils::install_material(context, bounce_material());
    constexpr real half = 128, thickness = 1;
    validation_extent_ = half + thickness + real(0.25);
    lite_walls(context, half, thickness);
    std::mt19937 rng(options.seed ^ options.body_count);
    std::uniform_real_distribution<float> size(0.125f, 0.5f);
    bodies_.clear();
    bodies_.reserve(options.body_count);
    std::vector<PlacedBody> placed;
    placed.reserve(options.body_count);
    for (u32 index = 0; index < options.body_count; ++index) {
        if (index % 2 == 0) {
            const real radius = size(rng);
            LiteCircleBodyDef definition;
            definition.position = non_overlapping_position(rng, half, radius, placed);
            definition.linear_velocity = velocity(rng, 400, 1000);
            track_body(context, scene_utils::lite_circle(context, definition, radius), radius);
            placed.push_back({definition.position, radius});
        } else {
            const Vec2 extents{size(rng), size(rng)};
            const real radius = extents.length();
            AabbBoxBodyDef definition;
            definition.position = non_overlapping_position(rng, half, radius, placed);
            definition.linear_velocity = velocity(rng, 300, 1000);
            track_body(context, scene_utils::aabb_box(context, definition, extents), radius);
            placed.push_back({definition.position, radius});
        }
    }
}

void ConvexPolygonBounceScenario::setup(Context& context, const Options& options) {
    context.settings.gravity = {};
    context.settings.contact.enable_dynamic_restitution = false;
    scene_utils::install_material(context, bounce_material());
    constexpr real half = 512, thickness = 1;
    validation_extent_ = half + thickness + real(0.25);
    rigid_walls(context, half, thickness);
    std::mt19937 rng(options.seed ^ (options.body_count * 0x9e3779b9u));
    std::uniform_real_distribution<float> rotation(0, static_cast<float>(2 * math::PI));
    std::uniform_real_distribution<float> angular(-0.5f, 0.5f);
    bodies_.clear();
    bodies_.reserve(options.body_count);
    std::vector<PlacedBody> placed;
    placed.reserve(options.body_count);
    for (u32 index = 0; index < options.body_count; ++index) {
        const auto first = polygon(rng);
        real radius = radius_of(first);
        RigidBodyDef definition;
        definition.position = non_overlapping_position(rng, 480, radius, placed);
        definition.rotation = rotation(rng);
        definition.linear_velocity = velocity(rng, 100, 300);
        definition.angular_velocity = angular(rng);
        const BodyHandle body = scene_utils::rigid_polygon(context, definition, first);
        if (index % 4 == 0) {
            const auto second = polygon(rng);
            const Vec2 offset{radius * real(0.25), 0};
            const ShapeHandle shape =
                scene_utils::polygon_shape(context, second, 1, offset, rotation(rng));
            scene_utils::attach_shape(context, body, shape);
            radius = std::max(radius, radius_of(second) + offset.x);
        }
        track_body(context, body, radius);
        placed.push_back({definition.position, radius});
    }
}

// Known issue: body escaped.
void MixedShapeBounceScenario::setup(Context& context, const Options& options) {
    context.settings.gravity = {};
    context.settings.contact.enable_dynamic_restitution = false;
    scene_utils::install_material(context, bounce_material());
    constexpr real half = 1024, thickness = 1;
    validation_extent_ = half + thickness + real(0.25);
    rigid_walls(context, half, thickness);
    std::mt19937 rng(options.seed ^ (options.body_count * 0x85ebca6bu));
    std::uniform_real_distribution<float> size(0.5f, 5);
    std::uniform_real_distribution<float> rotation(0, static_cast<float>(2 * math::PI));
    bodies_.clear();
    bodies_.reserve(options.body_count);
    std::vector<PlacedBody> placed;
    placed.reserve(options.body_count);
    for (u32 index = 0; index < options.body_count; ++index) {
        const u32 family = index % 9;
        const Vec2 linear_velocity = velocity(rng, 400, 500);
        if (family < 4) {
            const Vec2 extents{size(rng), size(rng)};
            const real radius = extents.length();
            AabbBoxBodyDef definition;
            definition.position = non_overlapping_position(rng, half, radius, placed);
            definition.linear_velocity = linear_velocity;
            track_body(context, scene_utils::aabb_box(context, definition, extents), radius);
            placed.push_back({definition.position, radius});
        } else if (family < 8) {
            const real radius = size(rng);
            LiteCircleBodyDef definition;
            definition.position = non_overlapping_position(rng, half, radius, placed);
            definition.linear_velocity = linear_velocity;
            track_body(context, scene_utils::lite_circle(context, definition, radius), radius);
            placed.push_back({definition.position, radius});
        } else {
            const auto points = polygon(rng);
            const real radius = radius_of(points);
            RigidBodyDef definition;
            definition.position = non_overlapping_position(rng, half, radius, placed);
            definition.rotation = rotation(rng);
            definition.linear_velocity = linear_velocity;
            track_body(context, scene_utils::rigid_polygon(context, definition, points), radius);
            placed.push_back({definition.position, radius});
        }
    }
}

} // namespace sqp::benchmark
