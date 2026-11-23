// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#pragma once

#include <cmath>

#include "sqphysics2d/build_config.h"
#include "sqphysics2d/detail/types.h"

namespace sqp {

// Two-dimensional vector with arithmetic in scalar type T.
template <class T> struct Vec2T {
    T x = T(0), y = T(0);
    constexpr Vec2T() = default;
    constexpr Vec2T(T x_value, T y_value) : x(x_value), y(y_value) {}
    template <class U>
    explicit constexpr Vec2T(const Vec2T<U>& other)
        : x(static_cast<T>(other.x)), y(static_cast<T>(other.y)) {}
    static constexpr Vec2T zero() {
        return {};
    }
    static constexpr Vec2T one() {
        return {T(1), T(1)};
    }
    constexpr Vec2T operator+() const {
        return *this;
    }
    constexpr Vec2T operator-() const {
        return {-x, -y};
    }
    constexpr Vec2T operator+(const Vec2T& v) const {
        return {x + v.x, y + v.y};
    }
    constexpr Vec2T operator-(const Vec2T& v) const {
        return {x - v.x, y - v.y};
    }
    constexpr Vec2T operator*(T value) const {
        return {x * value, y * value};
    }
    constexpr Vec2T operator/(T value) const {
        return {x / value, y / value};
    }
    Vec2T& operator+=(const Vec2T& v) {
        x += v.x;
        y += v.y;
        return *this;
    }
    Vec2T& operator-=(const Vec2T& v) {
        x -= v.x;
        y -= v.y;
        return *this;
    }
    Vec2T& operator*=(T value) {
        x *= value;
        y *= value;
        return *this;
    }
    Vec2T& operator/=(T value) {
        x /= value;
        y /= value;
        return *this;
    }
    bool operator==(const Vec2T&) const = default;
    // Squared Euclidean length, avoiding a square root.
    T length_squared() const {
        return x * x + y * y;
    }
    // Euclidean length.
    T length() const {
        return static_cast<T>(std::sqrt(length_squared()));
    }
    bool isfinite() const {
        return std::isfinite(x) && std::isfinite(y);
    }
};

template <class T> constexpr Vec2T<T> operator*(T s, const Vec2T<T>& v) {
    return v * s;
}

using Vec2f = Vec2T<real32>;
using Vec2d = Vec2T<real64>;
using Vec2 = Vec2T<real>;

// Axis-aligned box represented by inclusive minimum and maximum corners.
struct Aabb {
    Vec2 min, max;
    // Returns an inverted box suitable as an empty accumulation value.
    static Aabb empty();
    // True when both corners are finite and ordered component-wise.
    bool is_valid() const;
    // Midpoint and half-size of the box, respectively.
    Vec2 center() const;
    Vec2 extent() const;
    // Expands each side by the given distance.
    void expand(real);
    // Grows the box to include a point.
    void include(const Vec2&);
    // Inclusive overlap test.
    static bool overlap(const Aabb&, const Aabb&);
    // Smallest box containing both inputs.
    static Aabb merge(const Aabb&, const Aabb&);
    // Bounding box covering this box before and after a translation.
    static Aabb swept(const Aabb&, const Vec2&);
    bool operator==(const Aabb&) const = default;
};

// namespace sqp
} // namespace sqp

namespace sqp::math {

// Angle conversion constants.
inline constexpr real PI = real(3.1415926535897932384626433832795);
inline constexpr real64 PI_REAL64 = real64(3.1415926535897932384626433832795);
inline constexpr real INV_PI = real(1) / PI;
inline constexpr real DEG_PER_RAD = real(180) / PI;
inline constexpr real RAD_PER_DEG = PI / real(180);

// Converts simulation time units to seconds in the requested numeric type.
template <class T> constexpr T tu_to_seconds(tick_t value) {
    return static_cast<T>(value) / static_cast<T>(config::TU_PER_SECOND);
}
// Converts non-positive seconds to zero; positive values are rounded to nearest unit.
template <class T> constexpr tick_t seconds_to_tu(T value) {
    return value > T(0)
               ? static_cast<tick_t>(value * static_cast<T>(config::TU_PER_SECOND) + T(0.5))
               : 0;
}
template <class T> inline constexpr real radian_to_degree(T value) {
    return static_cast<real>(value) * DEG_PER_RAD;
}
template <class T> inline constexpr real degree_to_radian(T value) {
    return static_cast<real>(value) * RAD_PER_DEG;
}

// Thin standard-math wrappers provided in the sqp::math namespace.
template <typename T> inline T sqrt(T v) {
    return std::sqrt(v);
}

// These wrappers currently use the platform standard library, so transcendental
// results need not be bitwise identical across platforms.
template <typename T> inline T sin(T v) {
    return std::sin(v);
}
template <typename T> inline T cos(T v) {
    return std::cos(v);
}
template <typename T, typename U> inline auto pow(T base, U exponent) {
    return std::pow(base, exponent);
}

template <typename T> inline T abs(T v) {
    return std::abs(v);
}
template <typename T> inline T floor(T v) {
    return std::floor(v);
}
template <typename T> inline T ceil(T v) {
    return std::ceil(v);
}
template <typename T> inline T round(T v) {
    return std::round(v);
}
template <typename T> inline T fmod(T numerator, T denominator) {
    return std::fmod(numerator, denominator);
}
template <typename T> inline bool isfinite(T v) {
    return std::isfinite(v);
}

// Clamps a value to the closed unit interval.
template <typename T> inline constexpr real clamp01(T v) {
    v = v > T(1) ? T(1) : v;
    return v < T(0) ? T(0) : v;
}
// Wraps into [min_value, max_value); a zero-width range returns min_value.
template <typename T> inline T wrap(T v, T min_value, T max_value) {
    const T range = max_value - min_value;
    if (range == T(0))
        return min_value;
    v = fmod(v - min_value, range);
    if (v < T(0))
        v += range;
    return v + min_value;
}
template <typename T> inline T wrap_radian(T v) {
    return wrap(v, T(-PI), T(PI));
}
template <typename T> inline T wrap_degree(T v) {
    return wrap(v, T(-180), T(180));
}

// Dot/cross products, perpendiculars, normalization and rotation helpers.
template <class T> constexpr T dot(const Vec2T<T>& a, const Vec2T<T>& b) {
    return a.x * b.x + a.y * b.y;
}
template <class T> constexpr T cross(const Vec2T<T>& a, const Vec2T<T>& b) {
    return a.x * b.y - a.y * b.x;
}
template <class T> constexpr Vec2T<T> cross(T s, const Vec2T<T>& v) {
    return {-s * v.y, s * v.x};
}
template <class T> constexpr Vec2T<T> cross(const Vec2T<T>& v, T s) {
    return {s * v.y, -s * v.x};
}
template <class T> constexpr Vec2T<T> perp_ccw(const Vec2T<T>& v) {
    return {-v.y, v.x};
}
// Return +X for a vector too short to normalize reliably.
template <class T> inline Vec2T<T> safe_normalized(const Vec2T<T>& v) {
    const T len2 = v.length_squared();
    return len2 <= T(config::EPSILON) * T(config::EPSILON) ? Vec2T<T>{T(1), T(0)}
                                                           : v / std::sqrt(len2);
}
// Return zero for a vector too short to normalize reliably.
template <class T> inline Vec2T<T> normalized_or_zero(const Vec2T<T>& v) {
    const T len2 = v.length_squared();
    return len2 <= T(config::EPSILON) * T(config::EPSILON) ? Vec2T<T>{} : v / std::sqrt(len2);
}
template <class T> inline Vec2T<T> rotate_vector(const Vec2T<T>& v, T c, T s) {
    return {c * v.x - s * v.y, s * v.x + c * v.y};
}
template <class T> inline Vec2T<T> rotate_vector(const Vec2T<T>& v, T angle) {
    return rotate_vector(v, static_cast<T>(std::cos(angle)), static_cast<T>(std::sin(angle)));
}
template <class T> inline Vec2T<T> outward_edge_normal_ccw(const Vec2T<T>& a, const Vec2T<T>& b) {
    const auto edge = b - a;
    return safe_normalized(Vec2T<T>{edge.y, -edge.x});
}

// Generic bitmask helpers; masks are converted to the storage type.
template <typename Storage, typename Mask> constexpr bool has_any_flags(Storage flags, Mask mask) {
    return (flags & static_cast<Storage>(mask)) != Storage(0);
}
template <typename Storage, typename Mask> constexpr bool has_all_flags(Storage flags, Mask mask) {
    const Storage raw_mask = static_cast<Storage>(mask);
    return (flags & raw_mask) == raw_mask;
}
template <typename Storage, typename Mask> constexpr void enable_flags(Storage& flags, Mask mask) {
    flags = static_cast<Storage>(flags | static_cast<Storage>(mask));
}
template <typename Storage, typename Mask> constexpr void disable_flags(Storage& flags, Mask mask) {
    flags = static_cast<Storage>(flags & ~static_cast<Storage>(mask));
}
template <typename Storage, typename Mask>
constexpr void set_flags(Storage& flags, Mask mask, bool enabled) {
    enabled ? enable_flags(flags, mask) : disable_flags(flags, mask);
}

// Floors value/cell_size to an integer, including for negative values.
inline i32 floor_div(real value, real cell_size) {
    return static_cast<i32>(floor(value / cell_size));
}

static_assert(sizeof(Vec2) == sizeof(real) * 2);
static_assert(math::dot(Vec2{2, 3}, Vec2{4, 5}) == real(23));
static_assert(math::degree_to_radian(180) == math::PI);
static_assert(math::tu_to_seconds<real>(config::TU_PER_SECOND) == real(1));

} // namespace sqp::math
