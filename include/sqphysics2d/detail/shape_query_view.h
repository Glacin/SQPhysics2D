// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#pragma once

#include <variant>

#include "sqphysics2d/detail/data_storage.h"

namespace sqp::detail {

// Query-local geometry copy plus a borrowed ShapePool pointer.
// Circle geometry is centered at zero; local_center carries its body-local offset.
struct ShapeQueryView {
    using Geometry = std::variant<std::monostate, ConvexPoly, Circle>;

    const ShapePool* pool = nullptr;
    ShapeHandle handle;
    Geometry value;
    Vec2 local_center;
    real bound_radius_about_centroid = 0;
    real inner_bound_radius = 0;

    ShapeQueryView() = default;
    ShapeQueryView(const ShapePool& source, ShapeHandle id) : pool(&source), handle(id) {
        ShapeProperties properties;
        if (!source.get_properties(id, properties))
            return;
        local_center = properties.center_offset;
        bound_radius_about_centroid = properties.bound_radius;
        inner_bound_radius = properties.inner_bound_radius;
        if (id.kind == ShapeKind::Circle) {
            Circle circle;
            if (source.get_circle(id, circle))
                value = Circle({}, circle.radius);
        } else if (id.kind == ShapeKind::Polygon) {
            PolygonVertices points;
            if (source.get_polygon(id, points))
                value = ConvexPoly(points);
        }
    }

    const Geometry& geometry() const {
        return value;
    }
    bool is_valid() const {
        return !std::holds_alternative<std::monostate>(value);
    }
    // Conservative radius around the body origin for angular sweeps.
    real sweep_radius() const {
        return is_valid() ? local_center.length() + bound_radius_about_centroid : real(0);
    }
    ConvexPoly polygon_relative_to_body() const {
        PolygonVertices vertices;
        if (pool->get_polygon(handle, vertices))
            return ConvexPoly(vertices);
        return {};
    }
};

} // namespace sqp::detail

namespace sqp {

class PreparedPolygon {
public:
    // A query view already owns the same reconstructed polygon geometry.
    const ConvexPoly& get(const detail::ShapeQueryView& view) const {
        return std::get<ConvexPoly>(view.geometry());
    }

    // Cache the first polygon requested; create a new instance for a different shape.
    const ConvexPoly& get(const ShapePool& pool, ShapeHandle shape) {
        if (!ready_) {
            PolygonVertices vertices;
            if (pool.get_polygon(shape, vertices))
                polygon_ = ConvexPoly(vertices);
            ready_ = true;
        }
        return polygon_;
    }

private:
    ConvexPoly polygon_;
    bool ready_ = false;
};

} // namespace sqp
