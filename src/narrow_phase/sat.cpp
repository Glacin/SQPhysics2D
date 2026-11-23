// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "sat.h"

#include <algorithm>
#include <cstddef>
#include <limits>

namespace sqp {
namespace {

struct WorldPolygon {
    PolygonVertices vertices;
    Vec2 center;
};

WorldPolygon transform_poly(const ConvexPoly& poly, const Vec2& position, real rotation) {
    WorldPolygon result;
    const real cosine = math::cos(rotation);
    const real sine = math::sin(rotation);
    result.center = position + math::rotate_vector(poly.center_offset(), cosine, sine);
    for (const Vec2& point : poly.local)
        result.vertices.push_back(result.center + math::rotate_vector(point, cosine, sine));
    return result;
}

Vec2 world_center(const ConvexPoly& poly, const Vec2& position, real rotation) {
    return position +
           math::rotate_vector(poly.center_offset(), math::cos(rotation), math::sin(rotation));
}

void project(const PolygonVertices& vertices, const Vec2& axis, real& out_min, real& out_max) {
    out_min = out_max = math::dot(vertices[0], axis);
    for (std::size_t i = 1; i < vertices.size(); ++i) {
        const real value = math::dot(vertices[i], axis);
        out_min = std::min(out_min, value);
        out_max = std::max(out_max, value);
    }
}

SatHit scan_axes(const PolygonVertices& world_a, const PolygonVertices& world_b,
                 const Vec2& center_a, const Vec2& center_b) {
    SatHit hit;
    real best_depth = std::numeric_limits<real>::infinity();
    Vec2 best_axis(1, 0);
    bool best_ref_is_a = true;

    const auto test = [&](const PolygonVertices& reference, const PolygonVertices& incident,
                          bool ref_is_a) {
        for (std::size_t i = 0; i < reference.size(); ++i) {
            const Vec2 axis =
                math::outward_edge_normal_ccw(reference[i], reference[(i + 1) % reference.size()]);
            real min_reference, max_reference, min_incident, max_incident;
            project(reference, axis, min_reference, max_reference);
            project(incident, axis, min_incident, max_incident);
            // Measure both separating translations for contained projections.
            const real depth = std::min(max_reference - min_incident, max_incident - min_reference);
            if (depth < real(0))
                return false;
            if (depth < best_depth) {
                best_depth = depth;
                best_axis = axis;
                best_ref_is_a = ref_is_a;
            }
        }
        return true;
    };

    if (!test(world_a, world_b, true) || !test(world_b, world_a, false))
        return hit;

    Vec2 normal = best_ref_is_a ? best_axis : -best_axis;
    const Vec2 a_to_b = center_b - center_a;
    if (a_to_b.length_squared() > config::EPSILON * config::EPSILON &&
        math::dot(a_to_b, normal) < real(0))
        normal = -normal;

    hit.hit = true;
    hit.normal = normal;
    hit.depth = best_depth;
    hit.ref_is_A = best_ref_is_a;
    return hit;
}

SatHit intersect_world(const PolygonVertices& world_a, const PolygonVertices& world_b,
                       const Vec2& center_a, const Vec2& center_b) {
    SatHit hit = scan_axes(world_a, world_b, center_a, center_b);
    if (!hit.hit)
        return hit;

    const PolygonVertices& reference = hit.ref_is_A ? world_a : world_b;
    const Vec2 reference_normal = hit.ref_is_A ? hit.normal : -hit.normal;
    real best_alignment = -std::numeric_limits<real>::infinity();
    for (std::size_t i = 0; i < reference.size(); ++i) {
        const Vec2 edge_normal =
            math::outward_edge_normal_ccw(reference[i], reference[(i + 1) % reference.size()]);
        const real alignment = math::dot(edge_normal, reference_normal);
        if (alignment > best_alignment) {
            best_alignment = alignment;
            hit.ref_edge_index = static_cast<i32>(i);
        }
    }
    return hit;
}

} // namespace

PolygonVertices transform_poly_to_world(const ConvexPoly& poly, const Vec2& position,
                                        real rotation) {
    return transform_poly(poly, position, rotation).vertices;
}

SatHit sat_depth_query(const ConvexPoly& a, const Vec2& position_a, real rotation_a,
                       const ConvexPoly& b, const Vec2& position_b, real rotation_b) {
    if (a.local.size() < 3 || b.local.size() < 3)
        return {};
    const WorldPolygon world_a = transform_poly(a, position_a, rotation_a);
    const WorldPolygon world_b = transform_poly(b, position_b, rotation_b);
    return scan_axes(world_a.vertices, world_b.vertices, world_a.center, world_b.center);
}

SatHit sat_intersect(const ConvexPoly& a, const Vec2& position_a, real rotation_a,
                     const ConvexPoly& b, const Vec2& position_b, real rotation_b,
                     const PolygonVertices& world_a, const PolygonVertices& world_b) {
    if (a.local.size() < 3 || b.local.size() < 3 || world_a.size() != a.local.size() ||
        world_b.size() != b.local.size())
        return {};
    return intersect_world(world_a, world_b, world_center(a, position_a, rotation_a),
                           world_center(b, position_b, rotation_b));
}

} // namespace sqp
