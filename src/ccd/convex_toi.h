// SPDX-FileCopyrightText: 2022 Erin Catto
// SPDX-License-Identifier: MIT

// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

// Might need to change the algorithm because it is not designed for general dynamic-dynamic

// Known issue: need more iterations to converge for double precise.

#pragma once

#include "sqphysics2d/detail/math.h"
#include "sqphysics2d/detail/shape_query_view.h"
#include "sqphysics2d/shape.h"

namespace sqp {

enum class ConvexToiState : u8 {
    Hit,
    Separated,
    InitialOverlap,
    Failed,
};

struct ConvexToiOutput {
    ConvexToiState state = ConvexToiState::Separated;
    real fraction = real(1);
    Vec2 normal = Vec2(1, 0);
    u32 distance_iterations = 0;
    u32 root_iterations = 0;
};

struct ConvexDistanceOutput {
    real separation = 0;
    Vec2 normal = Vec2(1, 0); // A toward B
    Vec2 point_a;
    Vec2 point_b;
    u16 feature_a = 0;
    u16 feature_b = 0;
    bool valid = false;
    bool overlapped = false;
};

// Closest surface witnesses used by CCD and TGS speculative constraints.
ConvexDistanceOutput convex_signed_distance(const detail::ShapeQueryView& shape_a, Vec2 position_a,
                                            real rotation_a, const detail::ShapeQueryView& shape_b,
                                            Vec2 position_b, real rotation_b);

// Rotational continuous collision query derived from Box2D's time-of-impact
// separation-axis iteration. Angles are interpolated as scalars and are never
// normalized, so sweeps larger than one revolution retain their authored path.
ConvexToiOutput convex_time_of_impact(const detail::ShapeQueryView& shape_a, Vec2 position_a0,
                                      real rotation_a0, Vec2 position_a1, real rotation_a1,
                                      const detail::ShapeQueryView& shape_b, Vec2 position_b0,
                                      real rotation_b0, Vec2 position_b1, real rotation_b1,
                                      real target_penetration_m, real min_fraction = real(0),
                                      real max_fraction = real(1));

// Pool-backed entry points.
ConvexDistanceOutput convex_signed_distance(const ShapePool& pool, ShapeHandle shape_a,
                                            Vec2 position_a, real rotation_a, ShapeHandle shape_b,
                                            Vec2 position_b, real rotation_b);

ConvexToiOutput convex_time_of_impact(const ShapePool& pool, ShapeHandle shape_a, Vec2 position_a0,
                                      real rotation_a0, Vec2 position_a1, real rotation_a1,
                                      ShapeHandle shape_b, Vec2 position_b0, real rotation_b0,
                                      Vec2 position_b1, real rotation_b1, real target_penetration_m,
                                      real min_fraction = real(0), real max_fraction = real(1));

} // namespace sqp
