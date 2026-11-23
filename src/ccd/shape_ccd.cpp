// SPDX-FileCopyrightText: 2022 Erin Catto
// SPDX-License-Identifier: MIT

// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "shape_ccd.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

#include "ccd/convex_toi.h"

#include "narrow_phase/sat.h"
#include "sqphysics2d/body.h"
#include "sqphysics2d/detail/math.h"
#include "sqphysics2d/detail/shape_query_view.h"
#include "sqphysics2d/shape.h"

namespace sqp {

namespace {

struct SupportFeature {
    Vec2 first;
    Vec2 second;
};

Vec2 initial_overlap_normal(Vec2 point, const Aabb& box) {
    const real left = point.x - box.min.x;
    const real right = box.max.x - point.x;
    const real bottom = point.y - box.min.y;
    const real top = box.max.y - point.y;
    const real minimum = std::min({left, right, bottom, top});
    if (minimum == left) {
        return {-1, 0};
    }
    if (minimum == right) {
        return {1, 0};
    }
    if (minimum == bottom) {
        return {0, -1};
    }
    return {0, 1};
}

SupportFeature shape_support_feature(const detail::ShapeQueryView& shape, Vec2 position,
                                     real rotation, Vec2 direction) {
    if (!shape.is_valid())
        return {position, position};
    const Vec2 axis = math::safe_normalized(direction);
    if (std::holds_alternative<Circle>(shape.geometry())) {
        const Vec2 center = position + math::rotate_vector(shape.local_center, rotation);
        const Vec2 point = center + axis * std::get<Circle>(shape.geometry()).radius;
        return {point, point};
    }

    const auto& vertices = std::get<ConvexPoly>(shape.geometry()).local;
    if (vertices.empty())
        return {position, position};
    const real cosine = math::cos(rotation), sine = math::sin(rotation);
    const auto world_point = [&](std::size_t index) {
        return position + math::rotate_vector(vertices[index] + shape.local_center, cosine, sine);
    };
    real maximum = -std::numeric_limits<real>::infinity();
    std::size_t best_index = 0;
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        const real projection = math::dot(world_point(i), axis);
        if (projection > maximum) {
            maximum = projection;
            best_index = i;
        }
    }

    for (std::size_t i = 0; i < vertices.size(); ++i) {
        const Vec2 first = world_point(i);
        const Vec2 second = world_point((i + 1) % vertices.size());
        // A tiny angular error in a TOI normal can make the projection across a long
        // wall edge exceed a fixed tolerance and incorrectly select a distant corner.
        const real support_tolerance = std::max(real(1e-4), (second - first).length() * real(1e-4));
        if (maximum - math::dot(first, axis) <= support_tolerance &&
            maximum - math::dot(second, axis) <= support_tolerance)
            return {first, second};
    }
    const Vec2 best = world_point(best_index);
    return {best, best};
}

Vec2 feature_point_at(const SupportFeature& feature, Vec2 tangent, real coordinate) {
    const real first = math::dot(feature.first, tangent);
    const real second = math::dot(feature.second, tangent);
    const real span = second - first;
    if (math::abs(span) <= config::EPSILON)
        return feature.first;
    const real fraction = std::clamp((coordinate - first) / span, real(0), real(1));
    return feature.first + (feature.second - feature.first) * fraction;
}

Vec2 circle_center(const detail::ShapeQueryView& shape, Vec2 position, real rotation) {
    return position + math::rotate_vector(shape.local_center, rotation);
}

struct SweptShape {
    const detail::ShapeQueryView& shape;
    Vec2 position0;
    real rotation0;
    Vec2 position1;
    real rotation1;
    PreparedPolygon& prepared;

    Vec2 position(real time) const {
        return position0 + (position1 - position0) * time;
    }
    real rotation(real time) const {
        return rotation0 + (rotation1 - rotation0) * time;
    }
    const ConvexPoly& polygon() const {
        return prepared.get(shape);
    }
};

struct ShapeSweep {
    SweptShape a;
    SweptShape b;

    bool overlaps(real time, Vec2* out_normal = nullptr) const {
        const Vec2 position_a = a.position(time);
        const Vec2 position_b = b.position(time);
        const real rotation_a = a.rotation(time);
        const real rotation_b = b.rotation(time);
        if (std::holds_alternative<ConvexPoly>(a.shape.geometry()) &&
            std::holds_alternative<ConvexPoly>(b.shape.geometry())) {
            const SatHit hit = sat_depth_query(a.polygon(), position_a, rotation_a, b.polygon(),
                                               position_b, rotation_b);
            if (hit.hit && out_normal != nullptr)
                *out_normal = hit.normal;
            return hit.hit;
        }
        if (std::holds_alternative<Circle>(a.shape.geometry()) &&
            std::holds_alternative<Circle>(b.shape.geometry())) {
            const Vec2 delta = circle_center(b.shape, position_b, rotation_b) -
                               circle_center(a.shape, position_a, rotation_a);
            const real radius = std::get<Circle>(a.shape.geometry()).radius +
                                std::get<Circle>(b.shape.geometry()).radius;
            if (delta.length_squared() > radius * radius)
                return false;
            if (out_normal != nullptr)
                *out_normal = math::safe_normalized(delta);
            return true;
        }
        return false;
    }
};

void reset_sweep_result(real& out_toi, Vec2& out_normal) {
    out_toi = real(1);
    out_normal = Vec2(1, 0);
}

bool circle_pair_depth(const detail::ShapeQueryView& shape_a, Vec2 position_a, real rotation_a,
                       const detail::ShapeQueryView& shape_b, Vec2 position_b, real rotation_b,
                       PreparedPolygon& prepared_a, PreparedPolygon& prepared_b, real& depth,
                       Vec2& normal_ba) {
    const bool a_is_circle = std::holds_alternative<Circle>(shape_a.geometry());
    const bool b_is_circle = std::holds_alternative<Circle>(shape_b.geometry());
    if (a_is_circle && b_is_circle) {
        const Vec2 center_a = circle_center(shape_a, position_a, rotation_a);
        const Vec2 center_b = circle_center(shape_b, position_b, rotation_b);
        const real radius_a = std::get<Circle>(shape_a.geometry()).radius;
        const real radius_b = std::get<Circle>(shape_b.geometry()).radius;
        if (radius_a <= real(0) || radius_b <= real(0))
            return false;
        const Vec2 delta = center_b - center_a;
        const real distance_squared = delta.length_squared();
        const real radius_sum = radius_a + radius_b;
        if (distance_squared > radius_sum * radius_sum)
            return false;
        const real distance = distance_squared > config::EPSILON * config::EPSILON
                                  ? math::sqrt(distance_squared)
                                  : real(0);
        depth = std::max(real(0), radius_sum - distance);
        normal_ba = distance > config::EPSILON ? -delta / distance : Vec2(-1, 0);
        return true;
    }

    const auto& circle = a_is_circle ? shape_a : shape_b;
    const auto& polygon = a_is_circle ? shape_b : shape_a;
    const Vec2 circle_position = a_is_circle ? position_a : position_b;
    const real circle_rotation = a_is_circle ? rotation_a : rotation_b;
    const Vec2 polygon_position = a_is_circle ? position_b : position_a;
    const real polygon_rotation = a_is_circle ? rotation_b : rotation_a;
    const real radius = std::get<Circle>(circle.geometry()).radius;
    if (radius <= real(0))
        return false;
    PreparedPolygon& prepared = a_is_circle ? prepared_b : prepared_a;
    const ConvexPoly& polygon_geometry = prepared.get(polygon);
    if (polygon_geometry.local.size() < 3)
        return false;
    const PolygonVertices vertices =
        transform_poly_to_world(polygon_geometry, polygon_position, polygon_rotation);
    const Vec2 center = circle_center(circle, circle_position, circle_rotation);
    real max_separation = -std::numeric_limits<real>::infinity();
    i32 best_edge = 0;
    Vec2 best_normal(1, 0);
    const i32 count = static_cast<i32>(vertices.size());
    for (i32 i = 0; i < count; ++i) {
        const Vec2 edge_a = vertices[i];
        const Vec2 edge_b = vertices[(i + 1) % count];
        const Vec2 edge_normal = math::outward_edge_normal_ccw(edge_a, edge_b);
        const real separation = math::dot(edge_normal, center - edge_a);
        if (separation > radius)
            return false;
        if (separation > max_separation) {
            max_separation = separation;
            best_edge = i;
            best_normal = edge_normal;
        }
    }
    const Vec2 first = vertices[best_edge];
    const Vec2 edge = vertices[(best_edge + 1) % count] - first;
    const real edge_length_squared = edge.length_squared();
    const real along =
        edge_length_squared <= config::EPSILON * config::EPSILON
            ? real(0)
            : std::clamp(math::dot(center - first, edge) / edge_length_squared, real(0), real(1));
    const Vec2 closest = first + edge * along;
    Vec2 polygon_to_circle;
    if (max_separation > real(0)) {
        const Vec2 delta = center - closest;
        const real distance_squared = delta.length_squared();
        if (distance_squared > radius * radius)
            return false;
        const real distance = distance_squared > config::EPSILON * config::EPSILON
                                  ? math::sqrt(distance_squared)
                                  : real(0);
        polygon_to_circle = distance > config::EPSILON ? delta / distance : best_normal;
        depth = radius - distance;
    } else {
        polygon_to_circle = best_normal;
        depth = radius - max_separation;
    }
    depth = std::max(real(0), depth);
    normal_ba = a_is_circle ? polygon_to_circle : -polygon_to_circle;
    return true;
}

bool solve_point_circle_toi(Vec2 start, Vec2 end, Vec2 center, real radius, real& out_toi) {
    out_toi = real(1);
    if (radius <= real(0)) {
        return false;
    }
    const Vec2 d0 = start - center;
    const Vec2 velocity = end - start;
    const real c = d0.length_squared() - radius * radius;
    if (c <= config::EPSILON * config::EPSILON) {
        out_toi = real(0);
        return true;
    }
    const real a = velocity.length_squared();
    if (a <= config::EPSILON * config::EPSILON) {
        return false;
    }
    const real half_b = math::dot(d0, velocity);
    if (half_b >= real(0)) {
        return false;
    }
    const real discriminant = half_b * half_b - a * c;
    if (discriminant < real(0)) {
        return false;
    }
    const real q = -half_b + math::sqrt(discriminant);
    if (q <= config::EPSILON) {
        return false;
    }
    const real time = c / q;
    if (time < real(0) || time > real(1)) {
        return false;
    }
    out_toi = std::clamp(time, real(0), real(1));
    return true;
}

bool point_inside_poly(Vec2 point, const PolygonVertices& vertices) {
    for (size_t i = 0; i < vertices.size(); ++i) {
        if (math::dot(
                math::outward_edge_normal_ccw(vertices[i], vertices[(i + 1) % vertices.size()]),
                point - vertices[i]) > config::EPSILON) {
            return false;
        }
    }
    return true;
}

bool closest_poly_boundary_point(Vec2 point, const PolygonVertices& vertices, Vec2& out_point,
                                 Vec2& out_normal) {
    real best_distance2 = std::numeric_limits<real>::max();
    bool found = false;
    for (size_t i = 0; i < vertices.size(); ++i) {
        const Vec2 a = vertices[i];
        const Vec2 b = vertices[(i + 1) % vertices.size()];
        const Vec2 edge = b - a;
        const real edge_length2 = edge.length_squared();
        if (edge_length2 <= config::EPSILON * config::EPSILON) {
            continue;
        }
        const real u = std::clamp(math::dot(point - a, edge) / edge_length2, real(0), real(1));
        const Vec2 candidate = a + edge * u;
        const real distance2 = (point - candidate).length_squared();
        if (distance2 < best_distance2) {
            best_distance2 = distance2;
            out_point = candidate;
            out_normal = math::outward_edge_normal_ccw(a, b);
            found = true;
        }
    }
    return found;
}
} // namespace

// Fast analytic path for specific situations.
bool sweep_circle_vs_circle_toi(Vec2 a0, Vec2 a1, real radius_a, Vec2 b0, Vec2 b1, real radius_b,
                                bool enabled, real target_penetration_m, real& out_toi,
                                Vec2& out_normal) {
    reset_sweep_result(out_toi, out_normal);
    if (!enabled) {
        return false;
    }
    const real radius = std::max(real(0), radius_a + radius_b - target_penetration_m);
    if (radius <= real(0)) {
        return false;
    }
    const Vec2 start = a0 - b0;
    const Vec2 velocity = (a1 - a0) - (b1 - b0);
    const real c = start.length_squared() - radius * radius;
    if (c <= config::EPSILON * config::EPSILON) {
        out_toi = real(0);
        out_normal = math::safe_normalized(start);
        return true;
    }
    const real a = velocity.length_squared();
    if (a <= config::EPSILON * config::EPSILON) {
        return false;
    }
    const real b = math::dot(start, velocity);
    const real discriminant = b * b - a * c;
    if (b >= real(0) || discriminant < real(0)) {
        return false;
    }
    const real time = (-b - math::sqrt(discriminant)) / a;
    if (time < real(0) || time > real(1)) {
        return false;
    }
    out_toi = time;
    out_normal = math::safe_normalized(start + velocity * time);
    return true;
}
bool sweep_point_vs_aabb(Vec2 point, Vec2 delta, const Aabb& box, real& out_toi, Vec2& out_normal) {
    reset_sweep_result(out_toi, out_normal);
    if (!box.is_valid()) {
        return false;
    }

    if (point.x >= box.min.x - config::EPSILON && point.x <= box.max.x + config::EPSILON &&
        point.y >= box.min.y - config::EPSILON && point.y <= box.max.y + config::EPSILON) {
        out_toi = real(0);
        out_normal = initial_overlap_normal(point, box);
        return true;
    }

    real entry = -std::numeric_limits<real>::infinity();
    real exit = std::numeric_limits<real>::infinity();
    Vec2 entry_normal(1, 0);

    const auto update_axis = [&](real position, real velocity, real lower, real upper,
                                 Vec2 lower_normal, Vec2 upper_normal, real& axis_entry,
                                 real& axis_exit, Vec2& axis_normal) {
        if (math::abs(velocity) <= config::EPSILON) {
            return position >= lower - config::EPSILON && position <= upper + config::EPSILON;
        }
        real t0 = (lower - position) / velocity;
        real t1 = (upper - position) / velocity;
        Vec2 normal = lower_normal;
        if (t0 > t1) {
            std::swap(t0, t1);
            normal = upper_normal;
        }
        if (t0 > axis_entry) {
            axis_entry = t0;
            axis_normal = normal;
        }
        axis_exit = std::min(axis_exit, t1);
        return axis_entry <= axis_exit + config::EPSILON;
    };

    if (!update_axis(point.x, delta.x, box.min.x, box.max.x, {-1, 0}, {1, 0}, entry, exit,
                     entry_normal) ||
        !update_axis(point.y, delta.y, box.min.y, box.max.y, {0, -1}, {0, 1}, entry, exit,
                     entry_normal) ||
        exit < -config::EPSILON || entry > real(1) + config::EPSILON) {
        return false;
    }

    out_toi = std::clamp(entry, real(0), real(1));
    out_normal = entry_normal;
    return true;
}
bool sweep_circle_vs_aabb_toi(Vec2 circle0, Vec2 circle1, real circle_radius, const Aabb& box0,
                              const Aabb& box1, bool enabled, real target_penetration_m,
                              real& out_toi, Vec2& out_normal) {
    reset_sweep_result(out_toi, out_normal);
    if (!enabled || circle_radius <= real(0) || !box0.is_valid() || !box1.is_valid()) {
        return false;
    }

    const real radius = std::max(real(0), circle_radius - target_penetration_m);
    if (radius <= real(0)) {
        return false;
    }
    Aabb expanded = box0;
    expanded.expand(radius);
    const Vec2 relative_delta = (circle1 - circle0) - (box1.center() - box0.center());
    if (!sweep_point_vs_aabb(circle0, relative_delta, expanded, out_toi, out_normal))
        return false;
    const Vec2 center = circle0 + relative_delta * out_toi;
    const Vec2 corner{std::clamp(center.x, box0.min.x, box0.max.x),
                      std::clamp(center.y, box0.min.y, box0.max.y)};
    if ((center - corner).length_squared() <= radius * radius)
        return true;
    // The expanded box includes square corners outside the circle's rounded
    // Minkowski boundary. Refine such an entry against its nearest corner.
    if (!solve_point_circle_toi(circle0, circle0 + relative_delta, corner, radius, out_toi))
        return false;
    out_normal = math::safe_normalized(circle0 + relative_delta * out_toi - corner);
    return true;
}
bool sweep_aabb_vs_aabb_toi(const Aabb& a0, const Aabb& a1, const Aabb& b0, const Aabb& b1,
                            bool enabled, real target_penetration_m, real& out_toi,
                            Vec2& out_normal) {
    reset_sweep_result(out_toi, out_normal);
    if (!enabled || !a0.is_valid() || !a1.is_valid() || !b0.is_valid() || !b1.is_valid()) {
        return false;
    }

    const Vec2 extent_a = a0.extent();
    const real penetration = std::max(real(0), target_penetration_m);
    const Vec2 effective_extent(std::max(real(0), extent_a.x - penetration),
                                std::max(real(0), extent_a.y - penetration));
    Aabb expanded = b0;
    expanded.min -= effective_extent;
    expanded.max += effective_extent;
    const Vec2 relative_delta = (a1.center() - a0.center()) - (b1.center() - b0.center());
    return sweep_point_vs_aabb(a0.center(), relative_delta, expanded, out_toi, out_normal);
}
bool sweep_circle_vs_poly_toi(Vec2 circle0, Vec2 circle1, real circle_radius,
                              const ConvexPoly& poly, Vec2 poly_position, real poly_rotation0,
                              real poly_rotation1, Vec2 poly_linear_delta, bool enabled,
                              real target_penetration_m, real& out_toi, Vec2& out_normal,
                              Vec2& out_point) {
    reset_sweep_result(out_toi, out_normal);
    out_point = circle0;
    if (!enabled || circle_radius <= real(0) || poly.local.size() < 3) {
        return false;
    }

    const PolygonVertices vertices = transform_poly_to_world(poly, poly_position, poly_rotation0);
    const real angular_delta = math::abs(poly_rotation1 - poly_rotation0);
    real polygon_radius_squared = real(0);
    for (const Vec2& point : poly.local)
        polygon_radius_squared = std::max(polygon_radius_squared, point.length_squared());
    const real radius =
        std::max(real(0), circle_radius + angular_delta * math::sqrt(polygon_radius_squared) -
                              target_penetration_m);
    if (radius <= real(0)) {
        return false;
    }

    const Vec2 relative_end = circle1 - poly_linear_delta;
    const Vec2 velocity = relative_end - circle0;
    Vec2 initial_boundary;
    Vec2 initial_normal;
    if (closest_poly_boundary_point(circle0, vertices, initial_boundary, initial_normal)) {
        const Vec2 delta = circle0 - initial_boundary;
        if (point_inside_poly(circle0, vertices) ||
            delta.length_squared() <= radius * radius + config::EPSILON * config::EPSILON) {
            out_toi = real(0);
            out_normal = delta.length_squared() > config::EPSILON * config::EPSILON
                             ? math::safe_normalized(delta)
                             : initial_normal;
            out_point = initial_boundary;
            return true;
        }
    }
    if (velocity.length_squared() <= config::EPSILON * config::EPSILON) {
        return false;
    }

    bool found = false;
    real best_time = real(1);
    for (size_t i = 0; i < vertices.size(); ++i) {
        const Vec2 a = vertices[i];
        const Vec2 b = vertices[(i + 1) % vertices.size()];
        const Vec2 edge = b - a;
        const real edge_length2 = edge.length_squared();
        if (edge_length2 <= config::EPSILON * config::EPSILON) {
            continue;
        }
        const Vec2 normal = math::outward_edge_normal_ccw(a, b);
        const real separation0 = math::dot(normal, circle0 - a);
        const real separation1 = math::dot(normal, relative_end - a);
        const real denominator = separation1 - separation0;
        if (separation0 > radius && denominator < -config::EPSILON) {
            const real time = (radius - separation0) / denominator;
            if (time >= real(0) && time <= best_time) {
                const Vec2 center = circle0 + velocity * time;
                const real u = math::dot(center - a, edge) / edge_length2;
                if (u >= -config::EPSILON && u <= real(1) + config::EPSILON) {
                    found = true;
                    best_time = time;
                    out_normal = normal;
                    out_point = a + edge * std::clamp(u, real(0), real(1));
                }
            }
        }

        real vertex_time = real(1);
        if (solve_point_circle_toi(circle0, relative_end, a, radius, vertex_time) &&
            vertex_time < best_time) {
            found = true;
            best_time = vertex_time;
            const Vec2 center = circle0 + velocity * vertex_time;
            out_normal = math::safe_normalized(center - a);
            out_point = a;
        }
    }
    out_toi = best_time;
    return found;
}

CcdContactKinematics ccd_contact_kinematics(const detail::ShapeQueryView& shape_a, Vec2 position_a,
                                            real rotation_a, Vec2 linear_velocity_a,
                                            real angular_velocity_a,
                                            const detail::ShapeQueryView& shape_b, Vec2 position_b,
                                            real rotation_b, Vec2 linear_velocity_b,
                                            real angular_velocity_b, Vec2 normal_ba) {
    CcdContactKinematics result;
    result.normal_ab = -math::safe_normalized(normal_ba);
    const Vec2 tangent = math::perp_ccw(result.normal_ab);
    const SupportFeature feature_a =
        shape_support_feature(shape_a, position_a, rotation_a, result.normal_ab);
    const SupportFeature feature_b =
        shape_support_feature(shape_b, position_b, rotation_b, -result.normal_ab);
    const real a0 = math::dot(feature_a.first, tangent), a1 = math::dot(feature_a.second, tangent);
    const real b0 = math::dot(feature_b.first, tangent), b1 = math::dot(feature_b.second, tangent);
    const real a_min = std::min(a0, a1), a_max = std::max(a0, a1);
    const real b_min = std::min(b0, b1), b_max = std::max(b0, b1);
    const real overlap_min = std::max(a_min, b_min);
    const real overlap_max = std::min(a_max, b_max);

    std::array<real, 3> candidates{};
    if (overlap_min <= overlap_max) {
        candidates = {overlap_min, overlap_max, (overlap_min + overlap_max) * real(0.5)};
    } else if (a_max < b_min) {
        candidates = {a_max, b_min, (a_max + b_min) * real(0.5)};
    } else {
        candidates = {b_max, a_min, (b_max + a_min) * real(0.5)};
    }

    result.normal_speed = std::numeric_limits<real>::infinity();
    for (real coordinate : candidates) {
        const Vec2 point_a = feature_point_at(feature_a, tangent, coordinate);
        const Vec2 point_b = feature_point_at(feature_b, tangent, coordinate);
        const Vec2 velocity_a =
            linear_velocity_a + math::cross(angular_velocity_a, point_a - position_a);
        const Vec2 velocity_b =
            linear_velocity_b + math::cross(angular_velocity_b, point_b - position_b);
        const real speed = math::dot(velocity_b - velocity_a, result.normal_ab);
        if (speed < result.normal_speed) {
            result.point = (point_a + point_b) * real(0.5);
            result.normal_speed = speed;
        }
    }
    return result;
}

bool sweep_shape_pair_toi(const detail::ShapeQueryView& shape_a, Vec2 position_a0, real rotation_a0,
                          Vec2 position_a1, real rotation_a1, const detail::ShapeQueryView& shape_b,
                          Vec2 position_b0, real rotation_b0, Vec2 position_b1, real rotation_b1,
                          bool enabled, real target_penetration_m, real& out_toi, Vec2& out_normal,
                          u32* out_distance_iterations, u32* out_root_iterations,
                          bool* out_failed) {
    PreparedPolygon prepared_a, prepared_b;
    return sweep_shape_pair_toi(
        shape_a, position_a0, rotation_a0, position_a1, rotation_a1, shape_b, position_b0,
        rotation_b0, position_b1, rotation_b1, enabled, target_penetration_m, out_toi, out_normal,
        out_distance_iterations, out_root_iterations, out_failed, prepared_a, prepared_b);
}

bool sweep_shape_pair_toi(const detail::ShapeQueryView& shape_a, Vec2 position_a0, real rotation_a0,
                          Vec2 position_a1, real rotation_a1, const detail::ShapeQueryView& shape_b,
                          Vec2 position_b0, real rotation_b0, Vec2 position_b1, real rotation_b1,
                          bool enabled, real target_penetration_m, real& out_toi, Vec2& out_normal,
                          u32* out_distance_iterations, u32* out_root_iterations, bool* out_failed,
                          PreparedPolygon& prepared_a, PreparedPolygon& prepared_b) {
    reset_sweep_result(out_toi, out_normal);
    if (out_distance_iterations)
        *out_distance_iterations = 0;
    if (out_root_iterations)
        *out_root_iterations = 0;
    if (out_failed)
        *out_failed = false;
    if (!enabled) {
        return false;
    }
    if (!shape_a.is_valid() || !shape_b.is_valid()) {
        return false;
    }
    const ShapeSweep sweep{
        {shape_a, position_a0, rotation_a0, position_a1, rotation_a1, prepared_a},
        {shape_b, position_b0, rotation_b0, position_b1, rotation_b1, prepared_b}};
    const ConvexToiOutput result = convex_time_of_impact(
        shape_a, position_a0, rotation_a0, position_a1, rotation_a1, shape_b, position_b0,
        rotation_b0, position_b1, rotation_b1, target_penetration_m);
    if (out_distance_iterations)
        *out_distance_iterations = result.distance_iterations;
    if (out_root_iterations)
        *out_root_iterations = result.root_iterations;
    if (out_failed)
        *out_failed = result.state == ConvexToiState::Failed;
    if (result.state == ConvexToiState::Hit || result.state == ConvexToiState::InitialOverlap ||
        (result.state == ConvexToiState::Failed && result.fraction < real(1))) {
        out_toi = result.fraction;
        out_normal = result.normal;
        if (result.state == ConvexToiState::InitialOverlap) {
            Vec2 normal_ab;
            if (sweep.overlaps(real(0), &normal_ab))
                out_normal = -normal_ab;
        }
        if (result.state != ConvexToiState::InitialOverlap &&
            std::holds_alternative<ConvexPoly>(shape_a.geometry()) &&
            std::holds_alternative<ConvexPoly>(shape_b.geometry())) {
            Vec2 manifold_normal;
            const bool starts_overlapped = sweep.overlaps(real(0));
            if (!starts_overlapped) {
                real lower = real(0), upper = out_toi;
                bool bracketed = sweep.overlaps(upper, &manifold_normal);
                if (!bracketed) {
                    lower = out_toi;
                    const real motion_bound =
                        ((position_a1 - position_a0) - (position_b1 - position_b0)).length() +
                        math::abs(rotation_a1 - rotation_a0) * shape_a.sweep_radius() +
                        math::abs(rotation_b1 - rotation_b0) * shape_b.sweep_radius();
                    // A fast polygon can cross a thin wall between uniform samples.
                    // Probe densely next to the GJK estimate before widening the search.
                    real near_step = std::max(
                        config::EPSILON, real(0.002) / std::max(config::EPSILON, motion_bound));
                    while (near_step < real(1) / real(64) && !bracketed) {
                        upper = std::min(real(1), out_toi + near_step);
                        if (sweep.overlaps(upper, &manifold_normal))
                            bracketed = true;
                        else
                            lower = upper;
                        near_step *= real(2);
                    }
                    if (!bracketed)
                        lower = out_toi;
                    constexpr u32 samples = 64;
                    for (u32 sample = 1; sample <= samples && !bracketed; ++sample) {
                        upper = out_toi + (real(1) - out_toi) * real(sample) / real(samples);
                        if (sweep.overlaps(upper, &manifold_normal)) {
                            bracketed = true;
                            break;
                        }
                        lower = upper;
                    }
                }
                if (!bracketed)
                    return false;
                for (u32 iteration = 0; iteration < 24; ++iteration) {
                    const real middle = (lower + upper) * real(0.5);
                    if (sweep.overlaps(middle))
                        upper = middle;
                    else
                        lower = middle;
                }
                sweep.overlaps(upper, &manifold_normal);
                out_toi = upper;
                out_normal = -manifold_normal;
            }
        }
        return true;
    }
    return false;
}

bool sweep_shape_pair_reentry_toi(const detail::ShapeQueryView& shape_a, Vec2 position_a0,
                                  real rotation_a0, Vec2 position_a1, real rotation_a1,
                                  const detail::ShapeQueryView& shape_b, Vec2 position_b0,
                                  real rotation_b0, Vec2 position_b1, real rotation_b1,
                                  bool enabled, real target_penetration_m, real& out_toi,
                                  Vec2& out_normal) {
    PreparedPolygon prepared_a, prepared_b;
    return sweep_shape_pair_reentry_toi(shape_a, position_a0, rotation_a0, position_a1, rotation_a1,
                                        shape_b, position_b0, rotation_b0, position_b1, rotation_b1,
                                        enabled, target_penetration_m, out_toi, out_normal,
                                        prepared_a, prepared_b);
}

bool sweep_shape_pair_reentry_toi(const detail::ShapeQueryView& shape_a, Vec2 position_a0,
                                  real rotation_a0, Vec2 position_a1, real rotation_a1,
                                  const detail::ShapeQueryView& shape_b, Vec2 position_b0,
                                  real rotation_b0, Vec2 position_b1, real rotation_b1,
                                  bool enabled, real target_penetration_m, real& out_toi,
                                  Vec2& out_normal, PreparedPolygon& prepared_a,
                                  PreparedPolygon& prepared_b) {
    reset_sweep_result(out_toi, out_normal);
    if (!enabled || !shape_a.is_valid() || !shape_b.is_valid()) {
        return false;
    }
    const ShapeSweep sweep{
        {shape_a, position_a0, rotation_a0, position_a1, rotation_a1, prepared_a},
        {shape_b, position_b0, rotation_b0, position_b1, rotation_b1, prepared_b}};
    if (!sweep.overlaps(real(0)))
        return false;

    // Find an exit bracket with conservative advancement. This path is only
    // used for an initially overlapping, separating pair; it is deliberately
    // non-uniform and leaves the re-entry query to rotational TOI.
    const Vec2 relative_translation = (position_a1 - position_a0) - (position_b1 - position_b0);
    const real motion_bound = relative_translation.length() +
                              math::abs(rotation_a1 - rotation_a0) * shape_a.sweep_radius() +
                              math::abs(rotation_b1 - rotation_b0) * shape_b.sweep_radius();
    if (motion_bound <= config::EPSILON) {
        return false;
    }
    real lower = real(0);
    const real feature_distance =
        std::max(config::EPSILON,
                 std::min(shape_a.inner_bound_radius, shape_b.inner_bound_radius) * real(0.25));
    const real exit_advance = std::clamp(feature_distance / motion_bound, config::EPSILON, real(1));
    real upper = exit_advance;
    bool upper_overlaps = sweep.overlaps(upper);
    while (upper < real(1) && upper_overlaps) {
        lower = upper;
        upper = std::min(real(1), upper + exit_advance);
        upper_overlaps = sweep.overlaps(upper);
    }
    if (upper_overlaps) {
        return false;
    }
    for (u32 iteration = 0; iteration < 24; ++iteration) {
        const real middle = (lower + upper) * real(0.5);
        if (sweep.overlaps(middle)) {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    // Move beyond the numerical overlap band before starting the re-entry
    // TOI. At the binary exit itself GJK may intentionally classify a
    // touching simplex as InitialOverlap.
    upper = std::min(real(1), upper + real(16) * config::EPSILON / motion_bound);
    if (upper >= real(1) || sweep.overlaps(upper)) {
        return false;
    }
    // The overlap predicate changes at physical contact, while the TOI solver
    // accepts a small distance band around its target. Start outside that band
    // so the exit cannot be mistaken for a new impact.
    const real target_clearance =
        std::max(config::EPSILON * real(16),
                 std::min(shape_a.inner_bound_radius, shape_b.inner_bound_radius) * real(0.002));
    for (u32 retry = 0; retry < 64 && upper < real(1); ++retry) {
        if (convex_signed_distance(shape_a, sweep.a.position(upper), sweep.a.rotation(upper),
                                   shape_b, sweep.b.position(upper), sweep.b.rotation(upper))
                .separation > target_clearance)
            break;
        upper = std::min(real(1), upper + exit_advance);
        if (upper >= real(1) || sweep.overlaps(upper))
            return false;
    }
    if (target_penetration_m <= config::EPSILON) {
        const real feature_span =
            std::max(real(0.05),
                     std::min(shape_a.inner_bound_radius, shape_b.inner_bound_radius) * real(0.25));
        const u32 samples =
            std::clamp<u32>(static_cast<u32>(std::ceil(motion_bound / feature_span)), 64, 2048);
        real previous = upper;
        for (u32 sample = 1; sample <= samples; ++sample) {
            const real candidate =
                upper + (real(1) - upper) * real(sample) / real(samples);
            Vec2 manifold_normal;
            if (sweep.overlaps(candidate, &manifold_normal)) {
                real scan_lower = previous, hit = candidate;
                for (u32 iteration = 0; iteration < 24; ++iteration) {
                    const real middle = (scan_lower + hit) * real(0.5);
                    if (sweep.overlaps(middle))
                        hit = middle;
                    else
                        scan_lower = middle;
                }
                sweep.overlaps(hit, &manifold_normal);
                out_toi = hit;
                out_normal = -manifold_normal;
                return true;
            }
            previous = candidate;
        }
        return false;
    }
    ConvexToiOutput result = convex_time_of_impact(
        shape_a, position_a0, rotation_a0, position_a1, rotation_a1, shape_b, position_b0,
        rotation_b0, position_b1, rotation_b1, target_penetration_m, upper, real(1));
    for (u32 retry = 0; retry < 64 && result.state == ConvexToiState::Failed &&
                        result.fraction <= upper + config::EPSILON && upper < real(1);
         ++retry) {
        upper = std::min(real(1), upper + exit_advance * real(4));
        if (sweep.overlaps(upper))
            break;
        result = convex_time_of_impact(shape_a, position_a0, rotation_a0, position_a1, rotation_a1,
                                       shape_b, position_b0, rotation_b0, position_b1, rotation_b1,
                                       target_penetration_m, upper, real(1));
    }
    if (result.state == ConvexToiState::Hit ||
        (result.state == ConvexToiState::Failed && result.fraction < real(1))) {
        out_toi = result.fraction;
        out_normal = result.normal;
        return true;
    }
    return false;
}

bool sweep_shape_pair_overlap_approach_toi(const detail::ShapeQueryView& shape_a,
                                           const BodyState& start_a, const BodyState& end_a,
                                           const detail::ShapeQueryView& shape_b,
                                           const BodyState& start_b, const BodyState& end_b,
                                           real& out_toi, Vec2& out_normal) {
    reset_sweep_result(out_toi, out_normal);
    if (!shape_a.is_valid() || !shape_b.is_valid())
        return false;

    const real angular_a =
        math::abs(end_a.rotation - start_a.rotation) *
        (std::holds_alternative<ConvexPoly>(shape_a.geometry()) ? shape_a.sweep_radius()
                                                                : shape_a.local_center.length());
    const real angular_b =
        math::abs(end_b.rotation - start_b.rotation) *
        (std::holds_alternative<ConvexPoly>(shape_b.geometry()) ? shape_b.sweep_radius()
                                                                : shape_b.local_center.length());
    if (angular_a + angular_b <= config::EPSILON)
        return false;
    const Vec2 relative_translation =
        (end_a.position - start_a.position) - (end_b.position - start_b.position);
    const real motion_bound = relative_translation.length() + angular_a + angular_b;
    if (!std::isfinite(motion_bound) || motion_bound <= config::EPSILON)
        return false;

    // Contact-point velocity can report separation while a different rotating feature moves
    // into the obstacle. Track geometric depth from the shallowest pose instead.
    PreparedPolygon prepared_a, prepared_b;
    const ShapeSweep sweep{
        {shape_a, start_a.position, start_a.rotation, end_a.position, end_a.rotation, prepared_a},
        {shape_b, start_b.position, start_b.rotation, end_b.position, end_b.rotation, prepared_b}};
    const auto sample = [&](real fraction, real& depth, Vec2& normal) {
        const Vec2 position_a = sweep.a.position(fraction);
        const Vec2 position_b = sweep.b.position(fraction);
        const real rotation_a = sweep.a.rotation(fraction);
        const real rotation_b = sweep.b.rotation(fraction);

        if (std::holds_alternative<ConvexPoly>(shape_a.geometry()) &&
            std::holds_alternative<ConvexPoly>(shape_b.geometry())) {
            const SatHit hit = sat_depth_query(sweep.a.polygon(), position_a, rotation_a,
                                               sweep.b.polygon(), position_b, rotation_b);
            if (!hit.hit)
                return false;
            depth = std::max(real(0), hit.depth);
            normal = -hit.normal;
            return true;
        }
        return circle_pair_depth(shape_a, position_a, rotation_a, shape_b, position_b, rotation_b,
                                 prepared_a, prepared_b, depth, normal);
    };

    real minimum_depth = 0;
    Vec2 initial_normal;
    sample(real(0), minimum_depth, initial_normal);

    const real feature_span = std::max(
        real(0.05), std::min(shape_a.inner_bound_radius, shape_b.inner_bound_radius) * real(0.25));
    const u32 samples =
        static_cast<u32>(std::clamp(std::ceil(motion_bound / feature_span), real(8), real(2048)));
    constexpr real depth_tolerance = config::EPSILON * real(16);
    real previous = real(0);
    for (u32 index = 1; index <= samples; ++index) {
        const real fraction = real(index) / real(samples);
        real depth = 0;
        Vec2 normal;
        if (!sample(fraction, depth, normal)) {
            minimum_depth = real(0);
            previous = fraction;
            continue;
        }
        if (depth > minimum_depth + depth_tolerance) {
            real lower = previous, upper = fraction;
            for (u32 iteration = 0; iteration < 24; ++iteration) {
                const real middle = (lower + upper) * real(0.5);
                real middle_depth = 0;
                Vec2 middle_normal;
                if (sample(middle, middle_depth, middle_normal) &&
                    middle_depth > minimum_depth + depth_tolerance) {
                    upper = middle;
                    normal = middle_normal;
                } else {
                    lower = middle;
                }
            }
            out_toi = upper;
            out_normal = normal;
            return true;
        }
        minimum_depth = std::min(minimum_depth, depth);
        previous = fraction;
    }
    return false;
}

CcdContactKinematics ccd_contact_kinematics(const ShapePool& pool, ShapeHandle shape_a,
                                            Vec2 position_a, real rotation_a,
                                            Vec2 linear_velocity_a, real angular_velocity_a,
                                            ShapeHandle shape_b, Vec2 position_b, real rotation_b,
                                            Vec2 linear_velocity_b, real angular_velocity_b,
                                            Vec2 normal_ba) {
    return ccd_contact_kinematics(detail::ShapeQueryView(pool, shape_a), position_a, rotation_a,
                                  linear_velocity_a, angular_velocity_a,
                                  detail::ShapeQueryView(pool, shape_b), position_b, rotation_b,
                                  linear_velocity_b, angular_velocity_b, normal_ba);
}

bool sweep_shape_pair_toi(const ShapePool& pool, ShapeHandle shape_a, Vec2 position_a0,
                          real rotation_a0, Vec2 position_a1, real rotation_a1, ShapeHandle shape_b,
                          Vec2 position_b0, real rotation_b0, Vec2 position_b1, real rotation_b1,
                          bool enabled, real target_penetration_m, real& out_toi, Vec2& out_normal,
                          u32* out_distance_iterations, u32* out_root_iterations,
                          bool* out_failed) {
    return sweep_shape_pair_toi(detail::ShapeQueryView(pool, shape_a), position_a0, rotation_a0,
                                position_a1, rotation_a1, detail::ShapeQueryView(pool, shape_b),
                                position_b0, rotation_b0, position_b1, rotation_b1, enabled,
                                target_penetration_m, out_toi, out_normal, out_distance_iterations,
                                out_root_iterations, out_failed);
}

bool sweep_shape_pair_reentry_toi(const ShapePool& pool, ShapeHandle shape_a, Vec2 position_a0,
                                  real rotation_a0, Vec2 position_a1, real rotation_a1,
                                  ShapeHandle shape_b, Vec2 position_b0, real rotation_b0,
                                  Vec2 position_b1, real rotation_b1, bool enabled,
                                  real target_penetration_m, real& out_toi, Vec2& out_normal) {
    return sweep_shape_pair_reentry_toi(
        detail::ShapeQueryView(pool, shape_a), position_a0, rotation_a0, position_a1, rotation_a1,
        detail::ShapeQueryView(pool, shape_b), position_b0, rotation_b0, position_b1, rotation_b1,
        enabled, target_penetration_m, out_toi, out_normal);
}

bool sweep_shape_pair_overlap_approach_toi(const ShapePool& pool, ShapeHandle shape_a,
                                           const BodyState& start_a, const BodyState& end_a,
                                           ShapeHandle shape_b, const BodyState& start_b,
                                           const BodyState& end_b, real& out_toi,
                                           Vec2& out_normal) {
    return sweep_shape_pair_overlap_approach_toi(detail::ShapeQueryView(pool, shape_a), start_a,
                                                 end_a, detail::ShapeQueryView(pool, shape_b),
                                                 start_b, end_b, out_toi, out_normal);
}

} // namespace sqp
