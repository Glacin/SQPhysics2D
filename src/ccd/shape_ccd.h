// SPDX-FileCopyrightText: 2022 Erin Catto
// SPDX-License-Identifier: MIT

// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

// TODO: delete some function defines.

#pragma once

#include "sqphysics2d/body.h"
#include "sqphysics2d/detail/shape_query_view.h"
#include "sqphysics2d/shape.h"

namespace sqp {

// Sweep entry points use fixed EPSILON thresholds unless an epsilon is explicitly required.
// Motion arguments are grouped by shape, with start position/rotation then end position/rotation.

// Primitive CCD shortcuts. Normals point from shape B toward shape A.
// The caller must supply box0 and box1 for the same fixed-size, translating box.
// Reconstructed float bounds may differ slightly in extent, so this is not checked here.
bool sweep_circle_vs_aabb_toi(Vec2 circle0, Vec2 circle1, real circle_radius, const Aabb& box0,
                              const Aabb& box1, bool enabled, real target_penetration_m,
                              real& out_toi, Vec2& out_normal);

// The caller must supply each box's start/end bounds for pure translation at fixed size.
// Reconstructed float bounds may differ slightly in extent, so this is not checked here.
bool sweep_aabb_vs_aabb_toi(const Aabb& a0, const Aabb& a1, const Aabb& b0, const Aabb& b1,
                            bool enabled, real target_penetration_m, real& out_toi,
                            Vec2& out_normal);

bool sweep_circle_vs_circle_toi(Vec2 a0, Vec2 a1, real radius_a, Vec2 b0, Vec2 b1, real radius_b,
                                bool enabled, real target_penetration_m, real& out_toi,
                                Vec2& out_normal);

bool sweep_circle_vs_poly_toi(Vec2 circle0, Vec2 circle1, real circle_radius,
                              const ConvexPoly& poly, Vec2 poly_position, real poly_rotation0,
                              real poly_rotation1, Vec2 poly_linear_delta, bool enabled,
                              real target_penetration_m, real& out_toi, Vec2& out_normal,
                              Vec2& out_point);

struct CcdContactKinematics {
    Vec2 point;
    Vec2 normal_ab;
    real normal_speed = real(0);
};

CcdContactKinematics ccd_contact_kinematics(const detail::ShapeQueryView& shape_a, Vec2 position_a,
                                            real rotation_a, Vec2 linear_velocity_a,
                                            real angular_velocity_a,
                                            const detail::ShapeQueryView& shape_b, Vec2 position_b,
                                            real rotation_b, Vec2 linear_velocity_b,
                                            real angular_velocity_b, Vec2 normal_ba);

// The returned normal points from shape B toward shape A.
bool sweep_shape_pair_toi(const detail::ShapeQueryView& shape_a, Vec2 position_a0, real rotation_a0,
                          Vec2 position_a1, real rotation_a1, const detail::ShapeQueryView& shape_b,
                          Vec2 position_b0, real rotation_b0, Vec2 position_b1, real rotation_b1,
                          bool enabled, real target_penetration_m, real& out_toi, Vec2& out_normal,
                          u32* out_distance_iterations = nullptr,
                          u32* out_root_iterations = nullptr, bool* out_failed = nullptr);

bool sweep_shape_pair_toi(const detail::ShapeQueryView& shape_a, Vec2 position_a0, real rotation_a0,
                          Vec2 position_a1, real rotation_a1, const detail::ShapeQueryView& shape_b,
                          Vec2 position_b0, real rotation_b0, Vec2 position_b1, real rotation_b1,
                          bool enabled, real target_penetration_m, real& out_toi, Vec2& out_normal,
                          u32* out_distance_iterations, u32* out_root_iterations, bool* out_failed,
                          PreparedPolygon& prepared_a, PreparedPolygon& prepared_b);

// For a shapecast that starts overlapped but separating, find the first impact after the shapes
// have separated. The returned normal follows the same B-to-A contract.
bool sweep_shape_pair_reentry_toi(const detail::ShapeQueryView& shape_a, Vec2 position_a0,
                                  real rotation_a0, Vec2 position_a1, real rotation_a1,
                                  const detail::ShapeQueryView& shape_b, Vec2 position_b0,
                                  real rotation_b0, Vec2 position_b1, real rotation_b1,
                                  bool enabled, real target_penetration_m, real& out_toi,
                                  Vec2& out_normal);

bool sweep_shape_pair_reentry_toi(const detail::ShapeQueryView& shape_a, Vec2 position_a0,
                                  real rotation_a0, Vec2 position_a1, real rotation_a1,
                                  const detail::ShapeQueryView& shape_b, Vec2 position_b0,
                                  real rotation_b0, Vec2 position_b1, real rotation_b1,
                                  bool enabled, real target_penetration_m, real& out_toi,
                                  Vec2& out_normal, PreparedPolygon& prepared_a,
                                  PreparedPolygon& prepared_b);

// Find renewed penetration after an initial touch, including a short gap missed by re-entry.
// The normal points from B to A.
bool sweep_shape_pair_overlap_approach_toi(const detail::ShapeQueryView& shape_a,
                                           const BodyState& start_a, const BodyState& end_a,
                                           const detail::ShapeQueryView& shape_b,
                                           const BodyState& start_b, const BodyState& end_b,
                                           real& out_toi, Vec2& out_normal);

// Pool-backed entry points.
CcdContactKinematics ccd_contact_kinematics(const ShapePool& pool, ShapeHandle shape_a,
                                            Vec2 position_a, real rotation_a,
                                            Vec2 linear_velocity_a, real angular_velocity_a,
                                            ShapeHandle shape_b, Vec2 position_b, real rotation_b,
                                            Vec2 linear_velocity_b, real angular_velocity_b,
                                            Vec2 normal_ba);

bool sweep_shape_pair_toi(const ShapePool& pool, ShapeHandle shape_a, Vec2 position_a0,
                          real rotation_a0, Vec2 position_a1, real rotation_a1, ShapeHandle shape_b,
                          Vec2 position_b0, real rotation_b0, Vec2 position_b1, real rotation_b1,
                          bool enabled, real target_penetration_m, real& out_toi, Vec2& out_normal,
                          u32* out_distance_iterations = nullptr,
                          u32* out_root_iterations = nullptr, bool* out_failed = nullptr);

bool sweep_shape_pair_reentry_toi(const ShapePool& pool, ShapeHandle shape_a, Vec2 position_a0,
                                  real rotation_a0, Vec2 position_a1, real rotation_a1,
                                  ShapeHandle shape_b, Vec2 position_b0, real rotation_b0,
                                  Vec2 position_b1, real rotation_b1, bool enabled,
                                  real target_penetration_m, real& out_toi, Vec2& out_normal);

bool sweep_shape_pair_overlap_approach_toi(const ShapePool& pool, ShapeHandle shape_a,
                                           const BodyState& start_a, const BodyState& end_a,
                                           ShapeHandle shape_b, const BodyState& start_b,
                                           const BodyState& end_b, real& out_toi, Vec2& out_normal);

} // namespace sqp
