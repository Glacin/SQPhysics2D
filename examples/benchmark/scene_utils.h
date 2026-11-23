#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

#include "scenario_runner.h"

namespace sqp::benchmark::scene_utils {

inline PolygonVertices vertices(const std::vector<Vec2>& input) {
    PolygonVertices output;
    for (const Vec2 point : input)
        output.push_back(point);
    return output;
}

inline std::vector<Vec2> rectangle(real hx, real hy) {
    return {{-hx, -hy}, {hx, -hy}, {hx, hy}, {-hx, hy}};
}

inline void install_material(Context& context, const PhysicsMaterial& material) {
    if (!context.state.shape_pool.materials().set(0, material))
        throw std::runtime_error("material installation failed");
}

inline CommandResult command_result(Context& context, CommandId id) {
    CommandResult result;
    if (id == 0 || !context.world.step(0) || !context.world.take_command_result(id, result) ||
        result.status != CommandStatus::Applied)
        throw std::runtime_error("world setup command failed");
    return result;
}

inline void attach_shape(Context& context, BodyHandle body, ShapeHandle shape) {
    command_result(context, context.world.attach_shape(0, body, shape));
}

inline ShapeHandle polygon_shape(Context& context, const std::vector<Vec2>& points,
                                 real density = 1, Vec2 offset = {}, real rotation = 0) {
    PolygonShapeDef definition;
    definition.vertices = vertices(points);
    definition.local_offset = offset;
    definition.local_rotation = rotation;
    definition.density = density;
    definition.material_index = 0;
    const ShapeHandle shape = context.state.shape_pool.add_polygon(definition);
    if (!is_valid(shape))
        throw std::runtime_error("polygon shape creation failed");
    return shape;
}

inline ShapeHandle circle_shape(Context& context, real radius, real density = 1, Vec2 offset = {}) {
    CircleShapeDef definition;
    definition.radius = radius;
    definition.local_offset = offset;
    definition.density = density;
    definition.material_index = 0;
    const ShapeHandle shape = context.state.shape_pool.add_circle(definition);
    if (!is_valid(shape))
        throw std::runtime_error("circle shape creation failed");
    return shape;
}

inline BodyHandle rigid_polygon(Context& context, const RigidBodyDef& definition,
                                const std::vector<Vec2>& points, real density = 1, Vec2 offset = {},
                                real rotation = 0) {
    const ShapeHandle shape = polygon_shape(context, points, density, offset, rotation);
    const BodyHandle body = command_result(context, context.world.create_rigid(0, definition)).body;
    attach_shape(context, body, shape);
    return body;
}

inline BodyHandle lite_circle(Context& context, const LiteCircleBodyDef& definition, real radius,
                              real density = 1) {
    const ShapeHandle shape = circle_shape(context, radius, density);
    const BodyHandle body =
        command_result(context, context.world.create_lite_circle(0, definition, shape)).body;
    if (!is_valid(body))
        throw std::runtime_error("lite circle body creation failed");
    return body;
}

inline BodyHandle aabb_box(Context& context, const AabbBoxBodyDef& definition, Vec2 half_extents,
                           real density = 1) {
    if (density == 0) {
        RigidBodyDef fixed;
        fixed.position = definition.position;
        fixed.configs.register_type = BodyConfigs::RegisterType::ImmutableStatic;
        return rigid_polygon(context, fixed, rectangle(half_extents.x, half_extents.y), 0);
    }
    const ShapeHandle shape =
        polygon_shape(context, rectangle(half_extents.x, half_extents.y), density);
    const BodyHandle body =
        command_result(context, context.world.create_aabb_box(0, definition, shape)).body;
    if (!is_valid(body))
        throw std::runtime_error("AABB box body creation failed");
    return body;
}

inline bool finite(const BodyState& state) {
    return std::isfinite(state.position.x) && std::isfinite(state.position.y) &&
           std::isfinite(state.rotation) && std::isfinite(state.linear_velocity.x) &&
           std::isfinite(state.linear_velocity.y) && std::isfinite(state.angular_velocity);
}

inline Aabb body_bounds(const Context& context, BodyHandle body) {
    BodyState state;
    BodyProperties properties;
    if (!context.state.body_pool.get_state(body, state) ||
        !context.state.body_pool.get_properties(body, properties))
        return Aabb::empty();
    if (body.kind() == BodyConfigs::Kind::LiteCircle)
        return {state.position - Vec2(properties.radius, properties.radius),
                state.position + Vec2(properties.radius, properties.radius)};
    if (body.kind() == BodyConfigs::Kind::AabbBox)
        return {state.position - properties.half_extents, state.position + properties.half_extents};

    std::vector<ShapeHandle> shapes;
    if (!context.state.body_pool.get_attached_shapes(body, shapes))
        return Aabb::empty();
    const real cosine = math::cos(state.rotation), sine = math::sin(state.rotation);
    const auto transform = [&](Vec2 local) {
        local -= properties.local_center;
        return state.position + math::rotate_vector(local, cosine, sine);
    };
    Aabb bounds = Aabb::empty();
    for (const ShapeHandle shape : shapes) {
        if (shape.kind == ShapeKind::Circle) {
            Circle circle;
            if (!context.state.shape_pool.get_circle(shape, circle))
                continue;
            const Vec2 center = transform(circle.center);
            bounds.include(center - Vec2(circle.radius, circle.radius));
            bounds.include(center + Vec2(circle.radius, circle.radius));
        } else if (shape.kind == ShapeKind::Polygon) {
            PolygonVertices polygon;
            if (!context.state.shape_pool.get_polygon(shape, polygon))
                continue;
            for (const Vec2 point : polygon)
                bounds.include(transform(point));
        }
    }
    return bounds;
}

} // namespace sqp::benchmark::scene_utils
