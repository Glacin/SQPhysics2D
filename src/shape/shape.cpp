// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <cstddef>

#include "sqphysics2d/detail/math.h"
#include "sqphysics2d/shape.h"

namespace sqp {
namespace {

struct PolygonArea {
    real area = real(0);
    Vec2 centroid;
};

PolygonArea compute_polygon_area(const PolygonVertices& points) {
    PolygonArea result;
    if (points.size() < 3) {
        return result;
    }

    real twice_area = real(0);
    Vec2 weighted_centroid;
    for (std::size_t i = 0; i < points.size(); ++i) {
        const Vec2& a = points[i];
        const Vec2& b = points[(i + 1) % points.size()];
        const real edge_cross = math::cross(a, b);
        twice_area += edge_cross;
        weighted_centroid += (a + b) * edge_cross;
    }

    result.area = twice_area * real(0.5);
    if (result.area != real(0)) {
        result.centroid = weighted_centroid / (real(6) * result.area);
    }
    return result;
}

} // namespace

ConvexPoly::ConvexPoly(const PolygonVertices& points) {
    set_convex_points(points);
}

i32 ConvexPoly::vertex_count() const {
    return static_cast<i32>(local.size());
}

const Vec2& ConvexPoly::center_offset() const {
    return center_offset_;
}

void ConvexPoly::clear() {
    local.clear();
    center_offset_ = {};
}

bool ConvexPoly::set_convex_points(const PolygonVertices& points) {
    clear();
    if (points.overflowed()) {
        return false;
    }

    PolygonVertices cleaned = points;
    if (!sanitize_convex_ccw(cleaned)) {
        return false;
    }

    const PolygonArea area_result = compute_polygon_area(cleaned);
    if (area_result.area <= real(0)) {
        return false;
    }

    center_offset_ = area_result.centroid;
    for (const Vec2& point : cleaned) {
        local.push_back(point - center_offset_);
    }

    return true;
}

bool ConvexPoly::sanitize_convex_ccw(PolygonVertices& points) {
    if (points.overflowed() || points.size() < 3) {
        points.clear();
        return false;
    }

    const real eps2 = config::EPSILON * config::EPSILON;
    PolygonVertices cleaned;
    for (const Vec2& point : points) {
        if (!cleaned.empty() && (point - cleaned.back()).length_squared() <= eps2) {
            continue;
        }
        cleaned.push_back(point);
    }

    if (cleaned.size() >= 2 && (cleaned.front() - cleaned.back()).length_squared() <= eps2) {
        cleaned.pop_back();
    }
    if (cleaned.size() < 3) {
        points.clear();
        return false;
    }

    const PolygonArea area_result = compute_polygon_area(cleaned);
    if (math::abs(area_result.area) <= config::EPSILON || area_result.area < real(0)) {
        points.clear();
        return false;
    }

    bool removed = true;
    while (removed && cleaned.size() >= 3) {
        removed = false;
        const std::size_t count = cleaned.size();
        for (std::size_t i = 0; i < count; ++i) {
            const Vec2& a = cleaned[(i + count - 1) % count];
            const Vec2& b = cleaned[i];
            const Vec2& c = cleaned[(i + 1) % count];
            const real turn = math::cross(b - a, c - b);
            if (turn < -config::EPSILON) {
                points.clear();
                return false;
            }
            if (turn <= config::EPSILON || (c - b).length_squared() <= eps2) {
                cleaned.erase(i);
                removed = true;
                break;
            }
        }
    }

    if (cleaned.size() < 3) {
        points.clear();
        return false;
    }

    points = cleaned;
    return true;
}

Circle::Circle(const Vec2& center_value, real radius_value) {
    set(center_value, radius_value);
}

void Circle::clear() {
    center = {};
    radius = real(0);
}

void Circle::set(const Vec2& center_value, real radius_value) {
    center = center_value;
    radius = std::max(radius_value, real(0));
}

PolygonVertices::PolygonVertices(std::initializer_list<Vec2> points) {
    for (const auto& point : points)
        push_back(point);
}

void PolygonVertices::clear() {
    size_ = 0;
    overflowed_ = false;
}

void PolygonVertices::reserve(std::size_t) {}

void PolygonVertices::push_back(const Vec2& value) {
    if (size_ == values_.size()) {
        overflowed_ = true;
        return;
    }
    values_[size_++] = value;
}

void PolygonVertices::pop_back() {
    if (size_)
        --size_;
}

void PolygonVertices::erase(std::size_t index) {
    if (index >= size_)
        return;
    for (std::size_t i = index + 1; i < size_; ++i)
        values_[i - 1] = values_[i];
    --size_;
}

std::size_t PolygonVertices::size() const {
    return size_;
}

bool PolygonVertices::empty() const {
    return size_ == 0;
}

bool PolygonVertices::overflowed() const {
    return overflowed_;
}

Vec2* PolygonVertices::begin() {
    return values_.data();
}

Vec2* PolygonVertices::end() {
    return values_.data() + size_;
}

const Vec2* PolygonVertices::begin() const {
    return values_.data();
}

const Vec2* PolygonVertices::end() const {
    return values_.data() + size_;
}

Vec2& PolygonVertices::operator[](std::size_t i) {
    return values_[i];
}

const Vec2& PolygonVertices::operator[](std::size_t i) const {
    return values_[i];
}

Vec2& PolygonVertices::front() {
    return values_[0];
}

const Vec2& PolygonVertices::front() const {
    return values_[0];
}

Vec2& PolygonVertices::back() {
    return values_[size_ - 1];
}

const Vec2& PolygonVertices::back() const {
    return values_[size_ - 1];
}

} // namespace sqp
