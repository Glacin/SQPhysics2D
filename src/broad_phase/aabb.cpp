// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <limits>

#include "sqphysics2d/detail/math.h"

namespace sqp {

Aabb Aabb::empty() {
    const real inf = std::numeric_limits<real>::max();
    return {{inf, inf}, {-inf, -inf}};
}

bool Aabb::is_valid() const {
    return min.x <= max.x && min.y <= max.y;
}

Vec2 Aabb::center() const {
    return (min + max) * real(0.5);
}

Vec2 Aabb::extent() const {
    return (max - min) * real(0.5);
}

void Aabb::expand(real radius) {
    min -= Vec2{radius, radius};
    max += Vec2{radius, radius};
}

void Aabb::include(const Vec2& p) {
    min = {std::min(min.x, p.x), std::min(min.y, p.y)};
    max = {std::max(max.x, p.x), std::max(max.y, p.y)};
}

bool Aabb::overlap(const Aabb& a, const Aabb& b) {
    return a.min.x <= b.max.x && a.max.x >= b.min.x && a.min.y <= b.max.y && a.max.y >= b.min.y;
}

Aabb Aabb::merge(const Aabb& a, const Aabb& b) {
    return {{std::min(a.min.x, b.min.x), std::min(a.min.y, b.min.y)},
            {std::max(a.max.x, b.max.x), std::max(a.max.y, b.max.y)}};
}

Aabb Aabb::swept(const Aabb& a, const Vec2& delta) {
    Aabb result = a;
    delta.x > 0 ? result.max.x += delta.x : result.min.x += delta.x;
    delta.y > 0 ? result.max.y += delta.y : result.min.y += delta.y;
    return result;
}

} // namespace sqp
