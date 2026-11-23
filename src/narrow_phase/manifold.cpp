// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "manifold.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "narrow_phase/sat.h"

#include "sqphysics2d/detail/contact.h"
#include "sqphysics2d/detail/math.h"
#include "sqphysics2d/detail/shape_query_view.h"

namespace sqp {
namespace {

Vec2 shape_offset_world(const Vec2& position, real rotation, const Vec2& offset_from_body_com) {
    return position + math::rotate_vector(offset_from_body_com, rotation);
}

Vec2 closest_point_on_segment(const Vec2& point, const Vec2& a, const Vec2& b) {
    const Vec2 edge = b - a;
    const real length_squared = edge.length_squared();
    if (length_squared <= config::EPSILON * config::EPSILON) {
        return a;
    }
    const real t = std::clamp(math::dot(point - a, edge) / length_squared, real(0), real(1));
    return a + edge * t;
}

void set_circle_feature(ContactPoint& contact, bool is_a) {
    ContactFeatureRef& feature = is_a ? contact.key.feature_a : contact.key.feature_b;
    feature.type = ContactFeatureType::Circle;
    feature.index = 0;
    feature.auxiliary = 0;
}

i32 find_incident_edge_index(const PolygonVertices& polygon, const Vec2& reference_normal) {
    i32 best_index = 0;
    real best_dot = std::numeric_limits<real>::infinity();
    const i32 count = static_cast<i32>(polygon.size());
    for (i32 i = 0; i < count; ++i) {
        const Vec2 normal = math::outward_edge_normal_ccw(polygon[i], polygon[(i + 1) % count]);
        const real value = math::dot(normal, reference_normal);
        if (value < best_dot) {
            best_dot = value;
            best_index = i;
        }
    }
    return best_index;
}

struct ClipVertex {
    Vec2 point;
    ContactFeatureRef feature;
};

ContactFeatureRef make_vertex_feature(i32 index) {
    ContactFeatureRef feature;
    feature.type = ContactFeatureType::Vertex;
    feature.index = static_cast<u16>(index);
    return feature;
}

ContactFeatureRef make_edge_feature(i32 index, u16 aux = 0) {
    ContactFeatureRef feature;
    feature.type = ContactFeatureType::Edge;
    feature.index = static_cast<u16>(index);
    feature.auxiliary = aux;
    return feature;
}

i32 clip_segment_to_line(ClipVertex out[2], const ClipVertex in[2], const Vec2& normal, real offset,
                         real eps, i32 incident_edge, u16 clip_aux) {
    i32 out_count = 0;
    const real distance0 = math::dot(normal, in[0].point) - offset;
    const real distance1 = math::dot(normal, in[1].point) - offset;
    if (distance0 <= eps) {
        out[out_count++] = in[0];
    }
    if (distance1 <= eps) {
        out[out_count++] = in[1];
    }
    const bool crosses =
        (distance0 < -eps && distance1 > eps) || (distance0 > eps && distance1 < -eps);
    if (crosses) {
        const real t = distance0 / (distance0 - distance1);
        out[out_count++] = {in[0].point + (in[1].point - in[0].point) * t,
                            make_edge_feature(incident_edge, clip_aux)};
    }
    return out_count;
}

bool build_contacts_from_sat_pair(const SatHit& hit, const PolygonVertices& world_a,
                                  const PolygonVertices& world_b, ContactManifold& out) {
    const real normal_length_squared = hit.normal.length_squared();
    if (!hit.hit || !std::isfinite(hit.depth) || hit.depth < real(0) ||
        !std::isfinite(hit.normal.x) || !std::isfinite(hit.normal.y) ||
        !std::isfinite(normal_length_squared) || normal_length_squared <= real(0) ||
        hit.ref_edge_index < 0) {
        return false;
    }
    if (world_a.size() < 3 || world_b.size() < 3) {
        return false;
    }

    const Vec2 normal_ab = math::safe_normalized(hit.normal);
    const bool reference_is_a = hit.ref_is_A;
    const PolygonVertices& reference = reference_is_a ? world_a : world_b;
    const PolygonVertices& incident = reference_is_a ? world_b : world_a;
    const i32 reference_edge = hit.ref_edge_index;
    const i32 reference_count = static_cast<i32>(reference.size());
    if (reference_edge < 0 || reference_edge >= reference_count) {
        return false;
    }

    const Vec2 reference0 = reference[reference_edge];
    const Vec2 reference1 = reference[(reference_edge + 1) % reference_count];
    const Vec2 clip_normal = reference_is_a ? normal_ab : -normal_ab;
    const Vec2 tangent = math::safe_normalized(reference1 - reference0);
    const i32 incident_edge = find_incident_edge_index(incident, clip_normal);
    const i32 incident_count = static_cast<i32>(incident.size());
    const Vec2 incident0 = incident[incident_edge];
    const Vec2 incident1 = incident[(incident_edge + 1) % incident_count];

    const ClipVertex input[2] = {
        {incident0, make_vertex_feature(incident_edge)},
        {incident1, make_vertex_feature((incident_edge + 1) % incident_count)}};
    ClipVertex clipped1[2];
    ClipVertex clipped2[2];
    const Vec2 side_normal1 = -tangent;
    const real side_offset1 = math::dot(side_normal1, reference0);
    const i32 count1 = clip_segment_to_line(clipped1, input, side_normal1, side_offset1,
                                            config::EPSILON, incident_edge, 0);
    if (count1 <= 0) {
        return false;
    }
    const ClipVertex second_input[2] = {clipped1[0], count1 == 1 ? clipped1[0] : clipped1[1]};
    const Vec2 side_normal2 = tangent;
    const real side_offset2 = math::dot(side_normal2, reference1);
    const i32 count2 = clip_segment_to_line(clipped2, second_input, side_normal2, side_offset2,
                                            config::EPSILON, incident_edge, 1);
    if (count2 <= 0) {
        return false;
    }

    ClipVertex final_points[2];
    real final_separations[2];
    i32 final_count = 0;
    for (i32 i = 0; i < count2 && final_count < 2; ++i) {
        const real separation = math::dot(clip_normal, clipped2[i].point - reference0);
        // Keep both clipped endpoints, including a separated endpoint. Its
        // signed separation lets the temporal solver prevent future overlap.
        final_points[final_count] = clipped2[i];
        final_separations[final_count] = separation;
        ++final_count;
    }
    if (final_count == 0) {
        const Vec2 midpoint =
            ((reference0 + reference1) * real(0.5) + (incident0 + incident1) * real(0.5)) *
            real(0.5);
        final_points[0] = {midpoint, make_edge_feature(incident_edge, 0)};
        final_points[1] = {midpoint, make_edge_feature(incident_edge, 1)};
        final_separations[0] = math::dot(clip_normal, midpoint - reference0);
        final_separations[1] = final_separations[0];
        final_count = 2;
    } else if (final_count == 1) {
        final_points[1] = final_points[0];
        final_separations[1] = final_separations[0];
        final_count = 2;
    }

    out.count = 0;
    for (i32 i = 0; i < final_count && out.count < 2; ++i) {
        ContactPoint contact;
        // Store the midpoint between the incident point and its projection on
        // the reference plane. This is invariant under swapping A/B and gives
        // the solver symmetric lever arms for face contacts.
        contact.point = final_points[i].point - clip_normal * (final_separations[i] * real(0.5));

        contact.depth = -final_separations[i];

        if (reference_is_a) {
            contact.key.feature_a = make_edge_feature(reference_edge);
            contact.key.feature_b = final_points[i].feature;
        } else {
            contact.key.feature_b = make_edge_feature(reference_edge);
            contact.key.feature_a = final_points[i].feature;
        }
        out.points[out.count++] = contact;
    }
    out.normal = normal_ab;

    if (out.count == 2) {
        const Vec2 sort_axis = math::perp_ccw(normal_ab);
        if (math::dot(out.points[1].point, sort_axis) < math::dot(out.points[0].point, sort_axis))
            std::swap(out.points[0], out.points[1]);
    }
    return out.count > 0;
}

bool build_circle_circle_contacts(const Circle& a, const Vec2& center_a, const Circle& b,
                                  const Vec2& center_b, ContactManifold& out) {
    const real radius_a = a.radius;
    const real radius_b = b.radius;
    if (radius_a <= real(0) || radius_b <= real(0)) {
        return false;
    }
    const Vec2 delta = center_b - center_a;
    const real distance_squared = delta.length_squared();
    const real radius_sum = radius_a + radius_b;
    if (distance_squared > radius_sum * radius_sum) {
        return false;
    }
    const real distance = distance_squared > config::EPSILON * config::EPSILON
                              ? math::sqrt(distance_squared)
                              : real(0);
    const Vec2 normal = distance > config::EPSILON ? delta / distance : Vec2(1, 0);
    ContactPoint contact;
    contact.point = (center_a + normal * radius_a + center_b - normal * radius_b) * real(0.5);

    contact.depth = std::max(real(0), radius_sum - distance);
    set_circle_feature(contact, true);
    set_circle_feature(contact, false);
    out.points[0] = contact;
    out.count = 1;
    out.normal = normal;

    return true;
}

bool build_circle_poly_contacts(const Circle& circle, const Vec2& circle_center,
                                const ConvexPoly& polygon, const Vec2& polygon_position,
                                real polygon_rotation, bool circle_is_a, ContactManifold& out) {
    const real radius = circle.radius;
    if (radius <= real(0) || polygon.local.size() < 3) {
        return false;
    }
    const PolygonVertices world_polygon =
        transform_poly_to_world(polygon, polygon_position, polygon_rotation);
    const i32 vertex_count = static_cast<i32>(world_polygon.size());
    if (vertex_count < 3) {
        return false;
    }

    real max_separation = -std::numeric_limits<real>::infinity();
    i32 best_edge = 0;
    Vec2 best_normal(1, 0);
    for (i32 i = 0; i < vertex_count; ++i) {
        const Vec2 a = world_polygon[i];
        const Vec2 b = world_polygon[(i + 1) % vertex_count];
        const Vec2 normal = math::outward_edge_normal_ccw(a, b);
        const real separation = math::dot(normal, circle_center - a);
        if (separation > radius) {
            return false;
        }
        if (separation > max_separation) {
            max_separation = separation;
            best_edge = i;
            best_normal = normal;
        }
    }

    Vec2 closest = closest_point_on_segment(circle_center, world_polygon[best_edge],
                                            world_polygon[(best_edge + 1) % vertex_count]);
    Vec2 polygon_to_circle;
    real depth = real(0);
    const Vec2 center_delta = circle_center - closest;
    const real distance_squared = center_delta.length_squared();
    if (max_separation > real(0)) {
        if (distance_squared > radius * radius) {
            return false;
        }
        const real distance = distance_squared > config::EPSILON * config::EPSILON
                                  ? math::sqrt(distance_squared)
                                  : real(0);
        polygon_to_circle = distance > config::EPSILON ? center_delta / distance : best_normal;
        depth = radius - distance;
    } else {
        polygon_to_circle = best_normal;
        closest = circle_center - polygon_to_circle * max_separation;
        depth = radius - max_separation;
    }

    const Vec2 normal_ab = circle_is_a ? -polygon_to_circle : polygon_to_circle;
    ContactPoint contact;
    contact.point = closest;

    contact.depth = std::max(real(0), depth);
    if (circle_is_a) {
        set_circle_feature(contact, true);
        contact.key.feature_b.type = ContactFeatureType::Edge;
        contact.key.feature_b.index = static_cast<u16>(best_edge);
    } else {
        contact.key.feature_a.type = ContactFeatureType::Edge;
        contact.key.feature_a.index = static_cast<u16>(best_edge);
        set_circle_feature(contact, false);
    }
    out.points[0] = contact;
    out.count = 1;
    out.normal = normal_ab;

    return true;
}

} // namespace

bool build_contacts_from_shape_pair(const detail::ShapeQueryView& shape_a, const Vec2& position_a,
                                    real rotation_a, const detail::ShapeQueryView& shape_b,
                                    const Vec2& position_b, real rotation_b, ContactManifold& out) {
    out = ContactManifold();
    if (!shape_a.is_valid() || !shape_b.is_valid())
        return false;
    if (std::holds_alternative<ConvexPoly>(shape_a.geometry()) &&
        std::holds_alternative<ConvexPoly>(shape_b.geometry())) {
        if (std::get<ConvexPoly>(shape_a.geometry()).local.size() < 3 ||
            std::get<ConvexPoly>(shape_b.geometry()).local.size() < 3)
            return false;
        PreparedPolygon prepared_a, prepared_b;
        const ConvexPoly& polygon_a = prepared_a.get(shape_a);
        const ConvexPoly& polygon_b = prepared_b.get(shape_b);
        const auto world_a = transform_poly_to_world(polygon_a, position_a, rotation_a);
        const auto world_b = transform_poly_to_world(polygon_b, position_b, rotation_b);
        const SatHit hit = sat_intersect(polygon_a, position_a, rotation_a, polygon_b, position_b,
                                         rotation_b, world_a, world_b);
        return build_contacts_from_sat_pair(hit, world_a, world_b, out);
    }
    if (std::holds_alternative<Circle>(shape_a.geometry()) &&
        std::holds_alternative<Circle>(shape_b.geometry())) {
        const Vec2 center_a = shape_offset_world(position_a, rotation_a, shape_a.local_center);
        const Vec2 center_b = shape_offset_world(position_b, rotation_b, shape_b.local_center);
        return build_circle_circle_contacts(std::get<Circle>(shape_a.geometry()), center_a,
                                            std::get<Circle>(shape_b.geometry()), center_b, out);
    }
    if (std::holds_alternative<Circle>(shape_a.geometry())) {
        const ConvexPoly& polygon_b = std::get<ConvexPoly>(shape_b.geometry());
        if (polygon_b.local.size() < 3)
            return false;
        PreparedPolygon prepared_b;
        const Vec2 center_a = shape_offset_world(position_a, rotation_a, shape_a.local_center);
        return build_circle_poly_contacts(std::get<Circle>(shape_a.geometry()), center_a,
                                          prepared_b.get(shape_b), position_b,
                                          rotation_b, true, out);
    }
    const ConvexPoly& polygon_a = std::get<ConvexPoly>(shape_a.geometry());
    if (polygon_a.local.size() < 3)
        return false;
    PreparedPolygon prepared_a;
    const Vec2 center_b = shape_offset_world(position_b, rotation_b, shape_b.local_center);
    return build_circle_poly_contacts(std::get<Circle>(shape_b.geometry()), center_b,
                                      prepared_a.get(shape_a), position_a,
                                      rotation_a, false, out);
}

void ContactManifolds::clear() noexcept {
    keys.clear();
    points_x.clear();
    points_y.clear();
    depths.clear();
    counts.clear();
    normals_x.clear();
    normals_y.clear();
}

std::size_t ContactManifolds::size() const noexcept {
    return counts.size();
}

bool ContactManifolds::empty() const noexcept {
    return counts.empty();
}

bool ContactManifolds::append(const ContactManifold& manifold) {
    const std::size_t manifold_count = counts.size();
    const std::size_t point_count = manifold_count * 2;
    if (manifold.count < 1 || manifold.count > 2 || keys.size() != point_count ||
        points_x.size() != point_count || points_y.size() != point_count ||
        depths.size() != point_count || normals_x.size() != manifold_count ||
        normals_y.size() != manifold_count)
        return false;

    const auto grow = [](auto& column, std::size_t needed) {
        if (column.capacity() >= needed)
            return;
        const std::size_t capacity = column.capacity();
        const std::size_t doubled =
            capacity > column.max_size() / 2 ? column.max_size() : capacity * 2;
        column.reserve(std::max(needed, std::max(std::size_t(2), doubled)));
    };
    grow(keys, point_count + 2);
    grow(points_x, point_count + 2);
    grow(points_y, point_count + 2);
    grow(depths, point_count + 2);
    grow(counts, manifold_count + 1);
    grow(normals_x, manifold_count + 1);
    grow(normals_y, manifold_count + 1);

    const auto append_point = [&](const ContactPoint& point) {
        keys.push_back(point.key);
        points_x.push_back(point.point.x);
        points_y.push_back(point.point.y);
        depths.push_back(point.depth);
    };
    append_point(manifold.points[0]);
    append_point(manifold.count == 2 ? manifold.points[1] : ContactPoint{});
    counts.push_back(static_cast<u8>(manifold.count));
    normals_x.push_back(manifold.normal.x);
    normals_y.push_back(manifold.normal.y);
    return true;
}

ContactManifold ContactManifolds::at(std::size_t index) const {
    if (index >= size())
        throw std::out_of_range("ContactManifolds index");
    ContactManifold manifold;
    manifold.count = counts.at(index);
    manifold.normal = {normals_x.at(index), normals_y.at(index)};
    for (std::size_t slot = 0; slot < 2; ++slot) {
        const std::size_t point_index = index * 2 + slot;
        manifold.points[slot].key = keys.at(point_index);
        manifold.points[slot].point = {points_x.at(point_index), points_y.at(point_index)};
        manifold.points[slot].depth = depths.at(point_index);
    }
    return manifold;
}

} // namespace sqp
