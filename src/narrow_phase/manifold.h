// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#pragma once

#include "sqphysics2d/detail/contact.h"
#include "sqphysics2d/detail/shape_query_view.h"

namespace sqp {

// Builds the contact manifold for the shape pair at the supplied body poses.
bool build_contacts_from_shape_pair(const detail::ShapeQueryView& shape_a, const Vec2& position_a,
                                    real rotation_a, const detail::ShapeQueryView& shape_b,
                                    const Vec2& position_b, real rotation_b,
                                    ContactManifold& out_manifold);

} // namespace sqp
