#include "sqphysics2d/detail/data_storage.h"

#include <algorithm>
#include <limits>

namespace sqp {
namespace {

bool ray_segment(Vec2 from, Vec2 delta, Vec2 a, Vec2 b, real& fraction, Vec2& normal) {
    const Vec2 edge = b - a;
    const real denominator = math::cross(delta, edge);
    if (math::abs(denominator) <= config::EPSILON)
        return false;
    const Vec2 offset = a - from;
    const real t = math::cross(offset, edge) / denominator;
    const real u = math::cross(offset, delta) / denominator;
    if (t < real(0) || t > real(1) || u < real(0) || u > real(1))
        return false;
    fraction = t;
    normal = math::outward_edge_normal_ccw(a, b);
    if (math::dot(normal, delta) > real(0))
        normal = -normal;
    return true;
}

bool ray_circle(Vec2 from, Vec2 delta, Vec2 center, real radius, real& fraction, Vec2& normal) {
    const Vec2 offset = from - center;
    const real a = math::dot(delta, delta);
    if (a <= config::EPSILON * config::EPSILON)
        return false;
    const real b = real(2) * math::dot(offset, delta);
    const real c = math::dot(offset, offset) - radius * radius;
    const real discriminant = b * b - real(4) * a * c;
    if (discriminant < real(0))
        return false;
    const real root = math::sqrt(discriminant);
    const real t0 = (-b - root) / (real(2) * a);
    const real t1 = (-b + root) / (real(2) * a);
    const real t = t0 >= real(0) && t0 <= real(1) ? t0 : t1;
    if (t < real(0) || t > real(1))
        return false;
    fraction = t;
    normal = math::safe_normalized(from + delta * t - center);
    return true;
}

} // namespace

Aabb ShapePool::world_aabb(ShapeHandle handle, const BodyState& state) const {
    ShapeProperties properties;
    if (!get_properties(handle, properties))
        return Aabb::empty();
    if (handle.kind == ShapeKind::Circle) {
        Circle circle;
        if (!get_circle(handle, circle))
            return Aabb::empty();
        const Vec2 center =
            state.position + math::rotate_vector(properties.center_offset, state.rotation);
        return {center - Vec2(circle.radius, circle.radius),
                center + Vec2(circle.radius, circle.radius)};
    }
    PolygonVertices polygon;
    if (!get_polygon(handle, polygon))
        return Aabb::empty();
    Aabb box = Aabb::empty();
    const real cosine = math::cos(state.rotation), sine = math::sin(state.rotation);
    for (Vec2 point : polygon)
        box.include(state.position + math::rotate_vector(point, cosine, sine));
    return box;
}

real ShapePool::sweep_radius(ShapeHandle handle) const {
    ShapeProperties properties;
    return get_properties(handle, properties)
               ? properties.center_offset.length() + properties.bound_radius
               : real(0);
}

Aabb ShapePool::swept_world_aabb(ShapeHandle handle, const BodyState& start,
                                 const BodyState& end) const {
    if (!is_valid(handle))
        return Aabb::empty();
    Aabb bounds = Aabb::merge(world_aabb(handle, start), world_aabb(handle, end));
    if (math::abs(end.rotation - start.rotation) <= config::EPSILON)
        return bounds;
    const real radius = sweep_radius(handle);
    const Aabb rotational_bounds{{std::min(start.position.x, end.position.x) - radius,
                                  std::min(start.position.y, end.position.y) - radius},
                                 {std::max(start.position.x, end.position.x) + radius,
                                  std::max(start.position.y, end.position.y) + radius}};
    return Aabb::merge(bounds, rotational_bounds);
}

bool ShapePool::ray_cast(ShapeHandle handle, const BodyState& state, Vec2 from, Vec2 delta,
                         real& fraction, Vec2& normal) const {
    ShapeProperties properties;
    if (!get_properties(handle, properties))
        return false;
    if (handle.kind == ShapeKind::Circle) {
        Circle circle;
        if (!get_circle(handle, circle))
            return false;
        const Vec2 center =
            state.position + math::rotate_vector(properties.center_offset, state.rotation);
        return ray_circle(from, delta, center, circle.radius, fraction, normal);
    }
    PolygonVertices polygon;
    if (!get_polygon(handle, polygon))
        return false;
    bool hit = false;
    real best = std::numeric_limits<real>::max();
    Vec2 best_normal;
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const auto world_vertex = [&](std::size_t vertex) {
            return state.position + math::rotate_vector(polygon[vertex], state.rotation);
        };
        real candidate = real(0);
        Vec2 candidate_normal;
        if (ray_segment(from, delta, world_vertex(i), world_vertex((i + 1) % polygon.size()),
                        candidate, candidate_normal) &&
            candidate < best) {
            hit = true;
            best = candidate;
            best_normal = candidate_normal;
        }
    }
    if (hit) {
        fraction = best;
        normal = best_normal;
    }
    return hit;
}

} // namespace sqp
