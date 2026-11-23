// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

// SAT by myself, need to be improved later.

#pragma once

#include "sqphysics2d/shape.h"

namespace sqp {

struct SatHit {
    bool hit = false;
    Vec2 normal = Vec2(1, 0);
    real depth = real(0);
    bool ref_is_A = true;

    // -1 denotes a result without clipping-ready reference edge information.
    i32 ref_edge_index = -1;
};

PolygonVertices transform_poly_to_world(const ConvexPoly& poly, const Vec2& position,
                                        real rotation);

// Returns depth and normal without finding a clipping reference edge.
SatHit sat_depth_query(const ConvexPoly& a, const Vec2& position_a, real rotation_a,
                       const ConvexPoly& b, const Vec2& position_b, real rotation_b);

// Reuses world-space vertices already prepared by the caller.
SatHit sat_intersect(const ConvexPoly& a, const Vec2& position_a, real rotation_a,
                     const ConvexPoly& b, const Vec2& position_b, real rotation_b,
                     const PolygonVertices& world_a, const PolygonVertices& world_b);

} // namespace sqp
