// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <vector>

#include "sqphysics2d/body.h"
#include "sqphysics2d/detail/math.h"
#include "sqphysics2d/shape.h"

namespace sqp {

// Shape-level contact identity plus body motion versions used to reject stale pairs.
struct Pair {
    BodyHandle body_a, body_b;
    ShapeHandle shape_a, shape_b;
    u64 motion_version_a, motion_version_b;
    bool operator==(const Pair&) const = default;
};

// Which kind of geometric feature produced a contact point.
enum class ContactFeatureType : u8 {
    Edge = 0,
    Vertex = 1,
    Circle = 2,
};

// Identifies a specific feature (edge/vertex/circle) on a shape.
struct ContactFeatureRef {
    ContactFeatureType type = ContactFeatureType::Edge;
    u16 index = 0, auxiliary = 0;
    bool operator==(const ContactFeatureRef&) const = default;
};

// Persistent identity of a contact point within a shape pair.
struct ContactPointKey {
    ContactFeatureRef feature_a, feature_b;
    bool operator==(const ContactPointKey&) const = default;
};

// A contact point in the manifold's coordinate frame, with penetration depth.
struct ContactPoint {
    Vec2 point;
    real depth = 0;
    ContactPointKey key;
};

// Up to two contact points sharing a normal directed from A toward B.
// World callback manifolds use world-space points; TOI solver copies may use pair-local points.
struct ContactManifold {
    ContactPoint points[2];
    i32 count = 0;
    Vec2 normal{1, 0};
};

// Columnar manifold storage with two point slots per pair, including disabled contacts.
struct ContactManifolds {
    // size() == counts.size() * 2
    std::vector<ContactPointKey> keys;
    std::vector<real> points_x, points_y;
    std::vector<real> depths;

    // size() == counts.size()
    std::vector<u8> counts;
    std::vector<real> normals_x, normals_y;

    void clear() noexcept;
    std::size_t size() const noexcept;
    bool empty() const noexcept;
    // Appends a one- or two-point manifold; leaves the columns unchanged on invalid input.
    bool append(const ContactManifold&);
    // Reconstructs a single manifold. Throws std::out_of_range for an invalid index.
    ContactManifold at(std::size_t) const;
};

// Impulses applied at one contact point during a solve.
struct ContactPointImpulseReport {
    real normal_impulse = 0, tangent_impulse = 0;
};

// Per-contact summary of the solver's work for AfterSolve callbacks.
struct ContactSolveReport {
    ContactPointImpulseReport points[2];
    bool restitution_applied = false;
    real normal_impulse = 0, tangent_impulse = 0;
    real total_impulse = 0;
};

} // namespace sqp
