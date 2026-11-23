// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <cstddef>
#include <initializer_list>
#include <limits>

#include "sqphysics2d/detail/math.h"

namespace sqp {

// Fixed-capacity vertex container for convex polygons. An attempted insertion
// beyond capacity sets overflowed() and leaves existing vertices unchanged.
class PolygonVertices {
public:
    PolygonVertices() = default;
    PolygonVertices(std::initializer_list<Vec2>);

    void clear();
    // No-op hint for std::vector-like API compatibility; storage is fixed-size.
    void reserve(std::size_t);
    // Append unless full; an overflow attempt is remembered until clear().
    void push_back(const Vec2&);
    void pop_back();
    // Removes the indexed vertex, shifting later vertices down.
    void erase(std::size_t);

    std::size_t size() const;
    // Maximum number of vertices that can be stored.
    static constexpr std::size_t capacity() {
        return config::MAX_POLYGON_VERTICES;
    }

    bool empty() const;
    // True if insertion attempted to exceed capacity().
    bool overflowed() const;

    Vec2* begin();
    Vec2* end();
    const Vec2* begin() const;
    const Vec2* end() const;

    Vec2& operator[](std::size_t);
    const Vec2& operator[](std::size_t) const;
    Vec2& front();
    const Vec2& front() const;
    Vec2& back();
    const Vec2& back() const;

private:
    std::array<Vec2, config::MAX_POLYGON_VERTICES> values_{};
    std::size_t size_ = 0;
    bool overflowed_ = false;
};

// Parameters for creating or updating a circle shape.
// The shape is expressed in body-local metres; local_rotation is ignored for circles.
struct CircleShapeDef {
    real radius = 0;
    // Center offset from the body origin. Rotation has no effect on a circle.
    Vec2 local_offset;
    real local_rotation = 0;
    // Mass per unit area, used for mass and inertia.
    real density = 1;
    // Scales impact magnitude in restitution calculation.
    real restitution_impact_scale = 1;
    // Index into the ShapePool's MaterialPool.
    u32 material_index = 0;
};

// Parameters for creating or updating a convex polygon shape.
// Vertices must form a valid counter-clockwise convex polygon after the local transform.
struct PolygonShapeDef {
    // Convex hull points in local space.
    PolygonVertices vertices;
    // Shape offset and rotation relative to the body.
    Vec2 local_offset;
    real local_rotation = 0;
    // Mass per unit area, used for mass and inertia.
    real density = 1;
    // Scales impact magnitude in restitution calculation.
    real restitution_impact_scale = 1;
    // Index into the ShapePool's MaterialPool.
    u32 material_index = 0;
};

// Validated counter-clockwise convex polygon, recentered on its centroid.
struct ConvexPoly {
    // Counter-clockwise vertices relative to center_offset().
    PolygonVertices local;
    ConvexPoly() = default;
    // Builds from raw points; check vertex_count() to see whether they were accepted.
    explicit ConvexPoly(const PolygonVertices&);
    void clear();
    // Removes degeneracies, enforces counter-clockwise winding and stores the points.
    // Returns false and clears the polygon when the input is invalid.
    bool set_convex_points(const PolygonVertices&);
    i32 vertex_count() const;
    // Offset of the polygon's centroid from the local origin.
    const Vec2& center_offset() const;

private:
    Vec2 center_offset_{};
    // Cleans the vertex list in place into a counter-clockwise convex polygon.
    static bool sanitize_convex_ccw(PolygonVertices&);
};

// Circle represented by local center and nonnegative radius.
struct Circle {
    Circle() = default;
    Circle(const Vec2&, real);
    void clear();
    // Negative radii are clamped to zero.
    void set(const Vec2&, real);
    Vec2 center;
    real radius = 0;
};

// Selects the geometry family stored behind a ShapeHandle.
enum class ShapeKind : u8 { None = 0, Circle, Polygon };

inline constexpr u32 INVALID_SHAPE_INDEX = std::numeric_limits<u32>::max();

// Generation-checked reference to a shape in a ShapePool.
// The generation guards against use-after-remove: a stale handle whose slot
// was reused will fail is_valid().
struct ShapeHandle {
    u32 index = INVALID_SHAPE_INDEX;
    // 0 is never a live generation.
    u32 generation = 0;
    ShapeKind kind = ShapeKind::None;
    bool operator==(const ShapeHandle&) const = default;
};
// Cheap structural check (does not consult a pool; use ShapePool::is_valid
// to check the handle is still alive).
inline constexpr bool is_valid(ShapeHandle h) {
    return h.index != INVALID_SHAPE_INDEX && h.generation != 0 && h.kind != ShapeKind::None;
}

// Geometry, mass and material data derived from the stored shape definition.
struct ShapeProperties {
    // Shape centroid in body-local space.
    Vec2 center_offset;
    real area = 0;
    // Rotational inertia divided by mass.
    real inertia_per_unit_mass = 0;
    // Radius of a circle enclosing the shape.
    real bound_radius = 0;
    // Radius of a circle centered on the shape centroid and fully inside the shape.
    real inner_bound_radius = 0;
    // Mass per unit area and multiplier for dynamic restitution.
    real density = 1, restitution_impact_scale = 1;
    // Index into the owning ShapePool's MaterialPool.
    u32 material_index = 0;
};

struct BodyState;

class ShapePool;

} // namespace sqp
