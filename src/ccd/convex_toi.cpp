// SPDX-FileCopyrightText: 2022 Erin Catto
// SPDX-License-Identifier: MIT

// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "convex_toi.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include "sqphysics2d/detail/math.h"
#include "sqphysics2d/detail/shape_query_view.h"
#include "sqphysics2d/shape.h"

namespace sqp {
namespace {

struct DistanceProxy {
    std::array<Vec2, config::MAX_POLYGON_VERTICES> points{};
    u32 count = 0;
    real radius = real(0);
};

struct SimplexVertex {
    Vec2 point_a;
    Vec2 point_b;
    Vec2 point;
    real weight = real(0);
    u32 index_a = 0;
    u32 index_b = 0;
};

struct Simplex {
    std::array<SimplexVertex, 3> vertices{};
    u32 count = 0;
};

struct DistanceOutput {
    real distance = std::numeric_limits<real>::infinity();
    Vec2 normal = Vec2(1, 0);
    u32 iterations = 0;
    bool overlapped = false;
    Vec2 point_a;
    Vec2 point_b;
    u32 index_a = 0;
    u32 index_b = 0;
};

Vec2 transform_point(Vec2 position, real cosine, real sine, Vec2 point) {
    return position + math::rotate_vector(point, cosine, sine);
}

Vec2 inverse_rotate(real cosine, real sine, Vec2 vector) {
    return {cosine * vector.x + sine * vector.y, -sine * vector.x + cosine * vector.y};
}

DistanceProxy make_proxy(const detail::ShapeQueryView& shape) {
    DistanceProxy proxy;
    if (!shape.is_valid())
        return proxy;
    if (std::holds_alternative<Circle>(shape.geometry())) {
        proxy.points[0] = shape.local_center;
        proxy.count = 1;
        proxy.radius = std::get<Circle>(shape.geometry()).radius;
        return proxy;
    }
    if (std::holds_alternative<ConvexPoly>(shape.geometry())) {
        proxy.count = static_cast<u32>(std::get<ConvexPoly>(shape.geometry()).local.size());
        for (u32 i = 0; i < proxy.count; ++i) {
            proxy.points[i] = std::get<ConvexPoly>(shape.geometry()).local[i] + shape.local_center;
        }
    }
    return proxy;
}

u32 support(const DistanceProxy& proxy, Vec2 direction) {
    u32 best = 0;
    real best_value = math::dot(proxy.points[0], direction);
    for (u32 i = 1; i < proxy.count; ++i) {
        const real value = math::dot(proxy.points[i], direction);
        if (value > best_value) {
            best = i;
            best_value = value;
        }
    }
    return best;
}

Vec2 closest_point(const Simplex& simplex) {
    if (simplex.count == 1) {
        return simplex.vertices[0].point;
    }
    if (simplex.count == 2) {
        return simplex.vertices[0].point * simplex.vertices[0].weight +
               simplex.vertices[1].point * simplex.vertices[1].weight;
    }
    return {};
}

void solve_two(Simplex& simplex) {
    SimplexVertex& a = simplex.vertices[0];
    SimplexVertex& b = simplex.vertices[1];
    const Vec2 ab = b.point - a.point;
    const real region_b = math::dot(-a.point, ab);
    if (region_b <= real(0)) {
        a.weight = real(1);
        simplex.count = 1;
        return;
    }
    const real region_a = math::dot(b.point, ab);
    if (region_a <= real(0)) {
        b.weight = real(1);
        simplex.vertices[0] = b;
        simplex.count = 1;
        return;
    }
    const real inverse = real(1) / (region_a + region_b);
    a.weight = region_a * inverse;
    b.weight = region_b * inverse;
    simplex.count = 2;
}

void solve_three(Simplex& simplex) {
    SimplexVertex a = simplex.vertices[0];
    SimplexVertex b = simplex.vertices[1];
    SimplexVertex c = simplex.vertices[2];
    const Vec2 ab = b.point - a.point;
    const Vec2 ac = c.point - a.point;
    const Vec2 bc = c.point - b.point;
    const real d_ab_b = math::dot(-a.point, ab);
    const real d_ab_a = math::dot(b.point, ab);
    const real d_ac_c = math::dot(-a.point, ac);
    const real d_ac_a = math::dot(c.point, ac);
    const real d_bc_c = math::dot(-b.point, bc);
    const real d_bc_b = math::dot(c.point, bc);
    const real area = math::cross(ab, ac);
    const real area_a = area * math::cross(b.point, c.point);
    const real area_b = area * math::cross(c.point, a.point);
    const real area_c = area * math::cross(a.point, b.point);

    if (d_ab_b <= real(0) && d_ac_c <= real(0)) {
        a.weight = real(1);
        simplex.vertices[0] = a;
        simplex.count = 1;
        return;
    }
    if (d_ab_a > real(0) && d_ab_b > real(0) && area_c <= real(0)) {
        const real inverse = real(1) / (d_ab_a + d_ab_b);
        a.weight = d_ab_a * inverse;
        b.weight = d_ab_b * inverse;
        simplex.vertices[0] = a;
        simplex.vertices[1] = b;
        simplex.count = 2;
        return;
    }
    if (d_ac_a > real(0) && d_ac_c > real(0) && area_b <= real(0)) {
        const real inverse = real(1) / (d_ac_a + d_ac_c);
        a.weight = d_ac_a * inverse;
        c.weight = d_ac_c * inverse;
        simplex.vertices[0] = a;
        simplex.vertices[1] = c;
        simplex.count = 2;
        return;
    }
    if (d_ab_a <= real(0) && d_bc_c <= real(0)) {
        b.weight = real(1);
        simplex.vertices[0] = b;
        simplex.count = 1;
        return;
    }
    if (d_ac_a <= real(0) && d_bc_b <= real(0)) {
        c.weight = real(1);
        simplex.vertices[0] = c;
        simplex.count = 1;
        return;
    }
    if (d_bc_b > real(0) && d_bc_c > real(0) && area_a <= real(0)) {
        const real inverse = real(1) / (d_bc_b + d_bc_c);
        b.weight = d_bc_b * inverse;
        c.weight = d_bc_c * inverse;
        simplex.vertices[0] = b;
        simplex.vertices[1] = c;
        simplex.count = 2;
        return;
    }
    const real inverse = real(1) / (area_a + area_b + area_c);
    a.weight = area_a * inverse;
    b.weight = area_b * inverse;
    c.weight = area_c * inverse;
    simplex.vertices[0] = a;
    simplex.vertices[1] = b;
    simplex.vertices[2] = c;
    simplex.count = 3;
}

DistanceOutput shape_distance(const DistanceProxy& proxy_a, Vec2 position_a, real cosine_a,
                              real sine_a, const DistanceProxy& proxy_b, Vec2 position_b,
                              real cosine_b, real sine_b) {
    DistanceOutput output;
    Simplex simplex;
    Vec2 direction = position_b - position_a;
    if (direction.length_squared() <= config::EPSILON * config::EPSILON) {
        direction = Vec2(1, 0);
    }
    std::array<u32, 3> saved_a{};
    std::array<u32, 3> saved_b{};
    for (u32 iteration = 0; iteration < 20; ++iteration) {
        output.iterations = iteration + 1;
        const u32 saved_count = simplex.count;
        for (u32 i = 0; i < saved_count; ++i) {
            saved_a[i] = simplex.vertices[i].index_a;
            saved_b[i] = simplex.vertices[i].index_b;
        }
        if (simplex.count == 2) {
            solve_two(simplex);
        } else if (simplex.count == 3) {
            solve_three(simplex);
        }
        if (simplex.count == 3) {
            output.distance = real(0);
            output.overlapped = true;
            return output;
        }
        if (simplex.count > 0) {
            direction = -closest_point(simplex);
            if (direction.length_squared() <= config::EPSILON * config::EPSILON) {
                output.distance = real(0);
                output.overlapped = true;
                return output;
            }
        }
        SimplexVertex vertex;
        vertex.index_a = support(proxy_a, inverse_rotate(cosine_a, sine_a, direction));
        vertex.point_a =
            transform_point(position_a, cosine_a, sine_a, proxy_a.points[vertex.index_a]);
        vertex.index_b = support(proxy_b, inverse_rotate(cosine_b, sine_b, -direction));
        vertex.point_b =
            transform_point(position_b, cosine_b, sine_b, proxy_b.points[vertex.index_b]);
        vertex.point = vertex.point_a - vertex.point_b;
        bool duplicate = false;
        for (u32 i = 0; i < saved_count; ++i) {
            duplicate |= vertex.index_a == saved_a[i] && vertex.index_b == saved_b[i];
        }
        if (duplicate) {
            break;
        }
        simplex.vertices[simplex.count++] = vertex;
    }
    if (simplex.count == 2) {
        solve_two(simplex);
    }
    if (simplex.count == 1)
        simplex.vertices[0].weight = real(1);
    const Vec2 separation = closest_point(simplex);
    output.distance = separation.length();
    output.normal = math::safe_normalized(separation); // B toward A
    for (u32 i = 0; i < simplex.count; ++i) {
        output.point_a += simplex.vertices[i].point_a * simplex.vertices[i].weight;
        output.point_b += simplex.vertices[i].point_b * simplex.vertices[i].weight;
    }
    if (simplex.count > 0) {
        output.index_a = simplex.vertices[0].index_a;
        output.index_b = simplex.vertices[0].index_b;
    }
    return output;
}

real sweep_radius(const DistanceProxy& proxy) {
    real result = proxy.radius;
    for (u32 i = 0; i < proxy.count; ++i) {
        result = std::max(result, proxy.points[i].length() + proxy.radius);
    }
    return result;
}

// The minimum projected gap on a fixed world-space separating axis. This is a
// lower bound on the core-shape distance, even when the closest features change.
real axis_separation(const DistanceProxy& proxy_a, Vec2 position_a, real cosine_a, real sine_a,
                     const DistanceProxy& proxy_b, Vec2 position_b, real cosine_b, real sine_b,
                     Vec2 axis) {
    const u32 index_a = support(proxy_a, inverse_rotate(cosine_a, sine_a, axis));
    const u32 index_b = support(proxy_b, inverse_rotate(cosine_b, sine_b, -axis));
    const Vec2 point_a = transform_point(position_a, cosine_a, sine_a, proxy_a.points[index_a]);
    const Vec2 point_b = transform_point(position_b, cosine_b, sine_b, proxy_b.points[index_b]);
    return math::dot(point_b - point_a, axis);
}

} // namespace

ConvexDistanceOutput convex_signed_distance(const detail::ShapeQueryView& shape_a, Vec2 position_a,
                                            real rotation_a, const detail::ShapeQueryView& shape_b,
                                            Vec2 position_b, real rotation_b) {
    ConvexDistanceOutput output;
    const DistanceProxy proxy_a = make_proxy(shape_a);
    const DistanceProxy proxy_b = make_proxy(shape_b);
    if (proxy_a.count == 0 || proxy_b.count == 0)
        return output;
    const real cosine_a = math::cos(rotation_a);
    const real sine_a = math::sin(rotation_a);
    const real cosine_b = math::cos(rotation_b);
    const real sine_b = math::sin(rotation_b);
    const DistanceOutput distance = shape_distance(proxy_a, position_a, cosine_a, sine_a, proxy_b,
                                                   position_b, cosine_b, sine_b);
    output.valid = std::isfinite(distance.distance);
    output.overlapped = distance.overlapped;
    output.feature_a = static_cast<u16>(distance.index_a);
    output.feature_b = static_cast<u16>(distance.index_b);
    if (!output.valid)
        return output;
    output.normal = -distance.normal;
    output.separation = distance.distance - proxy_a.radius - proxy_b.radius;
    output.point_a = distance.point_a + output.normal * proxy_a.radius;
    output.point_b = distance.point_b - output.normal * proxy_b.radius;
    return output;
}

ConvexToiOutput convex_time_of_impact(const detail::ShapeQueryView& shape_a, Vec2 position_a0,
                                      real rotation_a0, Vec2 position_a1, real rotation_a1,
                                      const detail::ShapeQueryView& shape_b, Vec2 position_b0,
                                      real rotation_b0, Vec2 position_b1, real rotation_b1,
                                      real target_penetration_m, real min_fraction,
                                      real max_fraction) {
    ConvexToiOutput output;
    const DistanceProxy proxy_a = make_proxy(shape_a);
    const DistanceProxy proxy_b = make_proxy(shape_b);
    if (proxy_a.count == 0 || proxy_b.count == 0 || min_fraction > max_fraction) {
        output.state = ConvexToiState::Failed;
        return output;
    }

    const real total_radius = proxy_a.radius + proxy_b.radius;
    const real target =
        std::max(config::EPSILON, total_radius - std::max(real(0), target_penetration_m));

    const real tolerance = std::max(config::EPSILON * real(0.25), target * real(0.001));
    const Vec2 relative_translation = (position_a1 - position_a0) - (position_b1 - position_b0);
    const real angular_delta_a = rotation_a1 - rotation_a0;
    const real angular_delta_b = rotation_b1 - rotation_b0;
    const real angular_bound = math::abs(angular_delta_a) * sweep_radius(proxy_a) +
                               math::abs(angular_delta_b) * sweep_radius(proxy_b);
    const bool fixed_a = angular_delta_a == real(0);
    const bool fixed_b = angular_delta_b == real(0);
    const real fixed_angle_a = rotation_a0 + angular_delta_a * min_fraction;
    const real fixed_angle_b = rotation_b0 + angular_delta_b * min_fraction;
    const real fixed_cosine_a = fixed_a ? math::cos(fixed_angle_a) : real(0);
    const real fixed_sine_a = fixed_a ? math::sin(fixed_angle_a) : real(0);
    const real fixed_cosine_b = fixed_b ? math::cos(fixed_angle_b) : real(0);
    const real fixed_sine_b = fixed_b ? math::sin(fixed_angle_b) : real(0);

    const auto distance_at = [&](real time) {
        const Vec2 position_a = position_a0 + (position_a1 - position_a0) * time;
        const real angle_a = rotation_a0 + angular_delta_a * time;
        const real cosine_a = fixed_a ? fixed_cosine_a : math::cos(angle_a);
        const real sine_a = fixed_a ? fixed_sine_a : math::sin(angle_a);
        const Vec2 position_b = position_b0 + (position_b1 - position_b0) * time;
        const real angle_b = rotation_b0 + angular_delta_b * time;
        const real cosine_b = fixed_b ? fixed_cosine_b : math::cos(angle_b);
        const real sine_b = fixed_b ? fixed_sine_b : math::sin(angle_b);
        return shape_distance(proxy_a, position_a, cosine_a, sine_a, proxy_b, position_b, cosine_b,
                              sine_b);
    };

    real fraction = min_fraction;
    DistanceOutput distance = distance_at(fraction);
    output.distance_iterations += distance.iterations;
    const real separation0 = distance.distance - total_radius;
    if (separation0 <= config::EPSILON) {
        output.state = ConvexToiState::InitialOverlap;
        output.fraction = fraction;
        output.normal = distance.normal;
        return output;
    }

    const auto separation_at = [&](real time, Vec2 axis) {
        const Vec2 position_a = position_a0 + (position_a1 - position_a0) * time;
        const real angle_a = rotation_a0 + angular_delta_a * time;
        const Vec2 position_b = position_b0 + (position_b1 - position_b0) * time;
        const real angle_b = rotation_b0 + angular_delta_b * time;
        return axis_separation(proxy_a, position_a,
                               fixed_a ? fixed_cosine_a : math::cos(angle_a),
                               fixed_a ? fixed_sine_a : math::sin(angle_a), proxy_b, position_b,
                               fixed_b ? fixed_cosine_b : math::cos(angle_b),
                               fixed_b ? fixed_sine_b : math::sin(angle_b), axis);
    };

    // Advance a fixed Box2D-style separating axis using a bound on its closing
    // speed. A positive projected gap certifies that the full shapes are apart.
    // Refresh GJK when that certificate ends; this also finds transient hits
    // during sweeps whose final configuration is separated.

    // Might need more iterations for double precision.
    for (u32 outer = 0; outer < 20; ++outer) {
        const real separation = distance.distance;
        if (separation <= target + tolerance) {
            output.state = ConvexToiState::Hit;
            output.fraction = fraction;
            output.normal = distance.normal;
            return output;
        }
        const Vec2 axis = -distance.normal; // A toward B
        const real closing =
            std::max(real(0), math::dot(relative_translation, axis)) + angular_bound;
        if (closing <= config::EPSILON) {
            output.state = ConvexToiState::Separated;
            return output;
        }
        real safe_time = fraction;
        real safe_separation = separation;
        // Two polygon support scans cost about as much as the short GJK path in
        // this solver. Keep the scalar path for casts with a single-point proxy.
        bool scalar_valid = proxy_a.count == 1 || proxy_b.count == 1;
        if (scalar_valid) {
            safe_separation = separation_at(fraction, axis);
            const real axis_error = std::max(config::EPSILON * real(8), separation * real(0.0001));
            scalar_valid = std::isfinite(safe_separation) &&
                           math::abs(safe_separation - separation) <= axis_error;
            if (!scalar_valid)
                safe_separation = separation;
        }

        real candidate = fraction;
        real candidate_separation = safe_separation;
        DistanceOutput next;
        bool queried_distance = false;
        for (u32 scalar = 0; scalar < 4; ++scalar) {
            candidate = safe_time + (safe_separation - target) / closing;
            // EPSILON is a distance tolerance, not a fraction of the sweep.
            if (candidate <= safe_time || !std::isfinite(candidate)) {
                if (scalar_valid && safe_time > fraction) {
                    candidate = safe_time;
                    next = distance_at(candidate);
                    output.distance_iterations += next.iterations;
                    queried_distance = true;
                    break;
                }
                output.state = ConvexToiState::Failed;
                output.fraction = safe_time;
                output.normal = distance.normal;
                return output;
            }
            if (candidate > max_fraction) {
                output.state = ConvexToiState::Separated;
                output.fraction = max_fraction;
                output.normal = distance.normal;
                return output;
            }
            candidate_separation = scalar_valid ? separation_at(candidate, axis) : target;
            if (scalar_valid && std::isfinite(candidate_separation) &&
                candidate_separation > target + tolerance) {
                safe_time = candidate;
                safe_separation = candidate_separation;
                continue;
            }
            next = distance_at(candidate);
            output.distance_iterations += next.iterations;
            queried_distance = true;
            break;
        }
        if (!queried_distance) {
            // The scalar iteration budget is a performance limit, not a proof
            // of failure. Refresh GJK at the certified safe time.
            if (safe_time <= fraction) {
                output.state = ConvexToiState::Failed;
                output.fraction = fraction;
                output.normal = distance.normal;
                return output;
            }
            candidate = safe_time;
            next = distance_at(candidate);
            output.distance_iterations += next.iterations;
        }
        if (next.distance <= target + tolerance) {
            // First polish the cheap projected gap. Its root is only tentative:
            // the true Euclidean distance is checked before returning a hit.
            if (scalar_valid && candidate_separation < target - tolerance * real(4) &&
                next.distance < target - tolerance * real(4) && safe_separation > target) {
                real lower = safe_time, upper = candidate;
                real lower_value = safe_separation - target;
                real upper_value = candidate_separation - target;
                for (u32 root = 0; root < 50; ++root) {
                    ++output.root_iterations;
                    real trial = (lower + upper) * real(0.5);
                    if ((root & 1u) != 0u &&
                        math::abs(upper_value - lower_value) > config::EPSILON) {
                        trial = (lower * upper_value - upper * lower_value) /
                                (upper_value - lower_value);
                    }
                    if (trial <= lower || trial >= upper)
                        break;
                    const real value = separation_at(trial, axis) - target;
                    if (!std::isfinite(value))
                        break;
                    if (value > real(0)) {
                        lower = trial;
                        lower_value = value;
                    } else {
                        upper = trial;
                        upper_value = value;
                    }
                    if (math::abs(value) <= tolerance || upper - lower <= config::EPSILON)
                        break;
                }
                const DistanceOutput polished = distance_at(upper);
                output.distance_iterations += polished.iterations;
                if (polished.distance <= target + tolerance) {
                    output.state = ConvexToiState::Hit;
                    output.fraction = upper;
                    output.normal = math::dot(polished.normal, distance.normal) > real(0.5)
                                        ? polished.normal
                                        : distance.normal;
                    return output;
                }
            }
            // A changed feature or a numerical edge can invalidate the scalar
            // root. Fall back to the original full-distance bracket search.
            real lower = fraction;
            real upper = candidate;
            real lower_value = separation - target;
            real upper_value = next.distance - target;
            DistanceOutput root_distance = next;
            for (u32 root = 0; root < 50; ++root) {
                ++output.root_iterations;
                real trial;
                if ((root & 1u) == 0u && math::abs(upper_value - lower_value) > config::EPSILON) {
                    trial =
                        (lower * upper_value - upper * lower_value) / (upper_value - lower_value);
                } else {
                    trial = (lower + upper) * real(0.5);
                }
                trial = std::clamp(trial, lower, upper);
                root_distance = distance_at(trial);
                output.distance_iterations += root_distance.iterations;
                const real value = root_distance.distance - target;
                if (math::abs(value) <= tolerance || upper - lower <= config::EPSILON) {
                    upper = trial;
                    break;
                }
                if (value > real(0)) {
                    lower = trial;
                    lower_value = value;
                } else {
                    upper = trial;
                    upper_value = value;
                }
            }
            output.state = ConvexToiState::Hit;
            output.fraction = upper;
            output.normal =
                scalar_valid && math::dot(root_distance.normal, distance.normal) <= real(0.5)
                    ? distance.normal
                    : root_distance.normal;
            return output;
        }
        fraction = candidate;
        distance = next;
    }
    output.state = ConvexToiState::Failed;
    output.fraction = fraction;
    output.normal = distance.normal;
    return output;
}

ConvexDistanceOutput convex_signed_distance(const ShapePool& pool, ShapeHandle shape_a,
                                            Vec2 position_a, real rotation_a, ShapeHandle shape_b,
                                            Vec2 position_b, real rotation_b) {
    return convex_signed_distance(detail::ShapeQueryView(pool, shape_a), position_a, rotation_a,
                                  detail::ShapeQueryView(pool, shape_b), position_b, rotation_b);
}

ConvexToiOutput convex_time_of_impact(const ShapePool& pool, ShapeHandle shape_a, Vec2 position_a0,
                                      real rotation_a0, Vec2 position_a1, real rotation_a1,
                                      ShapeHandle shape_b, Vec2 position_b0, real rotation_b0,
                                      Vec2 position_b1, real rotation_b1, real target_penetration_m,
                                      real min_fraction, real max_fraction) {
    return convex_time_of_impact(detail::ShapeQueryView(pool, shape_a), position_a0, rotation_a0,
                                 position_a1, rotation_a1, detail::ShapeQueryView(pool, shape_b),
                                 position_b0, rotation_b0, position_b1, rotation_b1,
                                 target_penetration_m, min_fraction, max_fraction);
}

} // namespace sqp
