// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#pragma once

#include <memory>
#include <optional>
#include <span>
#include <type_traits>
#include <vector>

#include "sqphysics2d/body.h"
#include "sqphysics2d/joint.h"
#include "sqphysics2d/material.h"
#include "sqphysics2d/shape.h"

namespace sqp {

namespace detail { class StateHash; }

// Stores materials indexed by shapes. Index 0 is the default material.
class MaterialPool final {
public:
    MaterialPool();
    // Appends a material and returns its index.
    u32 add(const PhysicsMaterial&);
    // Copies the indexed material to the output; returns false for an invalid index.
    bool get(u32, PhysicsMaterial&) const;
    // Overwrites the indexed material; returns false for an invalid index.
    bool set(u32, const PhysicsMaterial&);
    // True if the index refers to an existing material.
    bool contains(u32) const;
    // Number of materials stored.
    u32 size() const;
    // Pre-allocates space for the given number of materials.
    void reserve(u32);

private:
    friend class ShapePool;
    // Restores the pool to its initial state (called by ShapePool::clear()).
    void reset();
    std::vector<PhysicsMaterial> materials_;
};

// Owns all collision shapes and their materials.
// Copies own independent shape storage while preserving handle values.
class ShapePool {
public:
    void append_state_hash(detail::StateHash&) const;
    ShapePool();
    ~ShapePool();
    ShapePool(const ShapePool&);
    ShapePool& operator=(const ShapePool&);
    ShapePool(ShapePool&&) noexcept;
    ShapePool& operator=(ShapePool&&) noexcept;

    // Removes all shapes (and resets materials).
    void clear();
    // Pre-allocates storage for the given number of shapes of each kind.
    void reserve(u32 circle_capacity, u32 polygon_capacity);

    // Creates a shape. Returns an invalid handle on failure
    // (e.g. invalid radius or non-convex/overflowed polygon).
    ShapeHandle add_circle(const CircleShapeDef&);
    ShapeHandle add_polygon(const PolygonShapeDef&);

    // Query shape geometry. Return false if the handle is invalid or of the wrong kind.
    bool get_circle(ShapeHandle, Circle&) const;
    bool get_polygon(ShapeHandle, PolygonVertices&) const;
    bool get_properties(ShapeHandle, ShapeProperties&) const;
    // Geometry queries use the pool's current shape data.
    Aabb world_aabb(ShapeHandle, const BodyState&) const;
    Aabb swept_world_aabb(ShapeHandle, const BodyState& start, const BodyState& end) const;
    real sweep_radius(ShapeHandle) const;
    bool ray_cast(ShapeHandle, const BodyState&, Vec2 from, Vec2 delta, real& fraction,
                  Vec2& normal) const;

    // Change density or restitution impact scale in place. Bodies using the shape need
    // BodyPool::refresh_properties() to pick up mass changes.
    bool set_density(ShapeHandle, real);
    bool set_restitution_impact_scale(ShapeHandle, real);

    // Material assignment / lookup for a shape.
    bool set_material_index(ShapeHandle, u32);
    bool get_material_index(ShapeHandle, u32&) const;
    bool get_material(ShapeHandle, PhysicsMaterial&) const;

    // Access to the material pool backing this shape pool.
    MaterialPool& materials();
    const MaterialPool& materials() const;

    // Replace a shape's definition in place, keeping the same handle.
    bool update_circle(ShapeHandle, const CircleShapeDef&);
    bool update_polygon(ShapeHandle, const PolygonShapeDef&);
    // Destroys a shape; its handle becomes invalid.
    bool remove(ShapeHandle);
    // True if the handle refers to a live shape in this pool.
    bool is_valid(ShapeHandle) const;

    ShapeKind kind(ShapeHandle) const;

    u32 shape_count() const;

private:
    friend class WorldCore;
    friend class WorldQuery;

    // Returns the validated convex polygon for a handle, or nullptr.
    const ConvexPoly* get_convex_polygon(ShapeHandle) const;
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// Stored real-valued columns. Fields not present for a body kind return an empty span.
enum class BodyScalarColumn : u8 {
    PositionX,
    PositionY,
    Rotation,
    VelocityX,
    VelocityY,
    AngularVelocity,
    ForceX,
    ForceY,
    Torque,
    Mass,
    InverseMass,
    Inertia,
    InverseInertia,
    DragMultiplier,
    LocalCenterX,
    LocalCenterY,
    BoundRadius,
    Density,
    Radius,
    HalfExtentX,
    HalfExtentY
};

// Owns all bodies. Copies preserve handle values and storage layout.
class IslandPool;
class JointPool;
class BodyPool {
public:
    void append_state_hash(detail::StateHash&) const;
    BodyPool();
    ~BodyPool();
    BodyPool(const BodyPool&);
    BodyPool& operator=(const BodyPool&);
    BodyPool(BodyPool&&) noexcept;
    BodyPool& operator=(BodyPool&&) noexcept;

    // A linked body requires the paired IslandPool for removal or clearing.
    bool clear();
    bool clear(IslandPool&);
    bool clear(JointPool&);
    bool clear(IslandPool&, JointPool&);
    // Pre-allocates storage per body kind.
    void reserve(u32 rigid_capacity, u32 circle_capacity, u32 box_capacity);

    // Create bodies. Lite circle / AABB bodies are bound to a single shape at creation (looked up
    // in the supplied ShapePool).
    // Rigid bodies start with no shapes; use attach_shape().
    // Return an invalid handle on failure.
    BodyHandle add_rigid(const RigidBodyDef& = {});
    BodyHandle add_lite_circle(const LiteCircleBodyDef&, ShapeHandle, const ShapePool&);
    BodyHandle add_aabb_box(const AabbBoxBodyDef&, ShapeHandle, const ShapePool&);

    // Destroys a body; its handle becomes invalid.
    bool remove(BodyHandle);
    bool remove(BodyHandle, IslandPool&);
    bool remove(BodyHandle, JointPool&);
    bool remove(BodyHandle, IslandPool&, JointPool&);
    // True if the handle refers to a live body in this pool.
    bool is_valid(BodyHandle) const;

    BodyConfigs::Kind kind(BodyHandle) const;
    u32 body_count() const;
    u32 rigid_count() const;
    u32 lite_circle_count() const;
    u32 aabb_box_count() const;

    // Fills the vector with handles of all live bodies.
    void get_handles(std::vector<BodyHandle>&) const;
    // Maps a live handle to its current index within that kind's dense columns.
    // Reacquire indices and spans after add, remove, reserve, clear, or moving the pool.
    std::optional<u32> dense_index(BodyHandle) const;
    // Returns no value for an invalid or unassigned body.
    std::optional<IslandId> island_of(BodyHandle) const;
    // Direct column writes bypass lock-aware setters.
    std::span<real> scalar_column(BodyConfigs::Kind, BodyScalarColumn);
    std::span<const real> scalar_column(BodyConfigs::Kind, BodyScalarColumn) const;
    std::span<BodyConfigs> configs_column(BodyConfigs::Kind);
    std::span<const BodyConfigs> configs_column(BodyConfigs::Kind) const;
    // Scheduler-owned version, indexed by the body's current dense index.
    std::span<u64> motion_version_column(BodyConfigs::Kind);
    std::span<const u64> motion_version_column(BodyConfigs::Kind) const;

    // Read / write a body's full state.
    bool get_state(BodyHandle, BodyState&) const;
    bool set_state(BodyHandle, const BodyState&);

    // Read derived / configured properties.
    bool get_properties(BodyHandle, BodyProperties&) const;
    bool get_configuration(BodyHandle, BodyConfiguration&) const;
    bool get_event_flags(BodyHandle, BodyConfigs::EventFlags&) const;
    // Lite bodies have the implicit Dynamic registration type.
    bool get_register_type(BodyHandle, BodyConfigs::RegisterType&) const;
    bool get_constant_force(BodyHandle, Vec2&, real& torque) const;

    // Setters. Each returns false if the handle is invalid.
    bool set_pose(BodyHandle, Vec2 position, real rotation);
    bool set_velocity(BodyHandle, Vec2 linear, real angular);
    bool set_drag_multiplier(BodyHandle, real);
    bool set_constant_force(BodyHandle, Vec2, real torque = 0);
    bool clear_forces(BodyHandle);
    bool set_ccd_override(BodyHandle, BodyConfigs::CcdOverride);
    bool set_collision_filter(BodyHandle, CollisionFilter);
    bool set_event_flags(BodyHandle, BodyConfigs::EventFlags);
    bool set_override_flags(BodyHandle, BodyConfigs::OverrideFlags);

    // Attach / detach a shape to a rigid body. Mass properties are
    // recomputed from the attached shapes.
    bool attach_shape(BodyHandle, ShapeHandle, const ShapePool&);
    bool detach_shape(BodyHandle, ShapeHandle, const ShapePool&);

    // Recomputes mass/inertia/bounds, e.g. after changing a shape's density.
    bool refresh_properties(BodyHandle, const ShapePool&);
    // Lists the shapes attached to a body.
    bool get_attached_shapes(BodyHandle, std::vector<ShapeHandle>&) const;
    // Finds the body that owns a given shape.
    bool get_shape_owner(ShapeHandle, BodyHandle&) const;
    // Gets a shape's offset from the body's center of mass.
    bool shape_offset_from_com(BodyHandle, ShapeHandle, const ShapePool&, real&, real&) const;

private:
    friend class IslandPool;
    friend class JointPool;
    bool assign_island(BodyHandle, IslandId);
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

static_assert(std::is_copy_constructible_v<BodyPool>);
static_assert(std::is_copy_constructible_v<ShapePool>);

// Owns joints registered with one BodyPool. Remove its joints before destroying either pool.
// The default World pipeline currently stores but does not solve these joints.
class JointPool {
public:
    void append_state_hash(detail::StateHash&) const;
    JointPool();
    ~JointPool();
    JointPool(const JointPool&) = delete;
    JointPool& operator=(const JointPool&) = delete;
    JointPool(JointPool&&) noexcept;
    JointPool& operator=(JointPool&&) noexcept;

    // Clone with the body pool that owns the copied joint links.
    void copy_from(const JointPool&, BodyPool& destination_bodies);

    JointHandle create(BodyPool&, const PinJointDef&);
    JointHandle create(BodyPool&, const SliderJointDef&);
    JointHandle create(BodyPool&, const WeldJointDef&);
    JointHandle create(BodyPool&, const DistanceJointDef&);
    JointHandle create(BodyPool&, const SpringJointDef&);

    bool is_valid(JointHandle) const;
    // Returned pointers can be invalidated by creation or deletion in this pool.
    const JointRegistration* get(JointHandle) const;
    const PinJoint* get_pin(JointHandle) const;
    const SliderJoint* get_slider(JointHandle) const;
    const WeldJoint* get_weld(JointHandle) const;
    const DistanceJoint* get_distance(JointHandle) const;
    const SpringJoint* get_spring(JointHandle) const;

    // Traversal is invalidated by any mutation to this body's joint chain.
    JointHandle first_attached(const BodyPool&, BodyHandle) const;
    JointHandle next_attached(BodyHandle, JointHandle) const;
    bool destroy(BodyPool&, JointHandle);
    bool destroy_attached(BodyPool&, BodyHandle);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// Groups movable bodies that advance together and share a scheduler time.
class IslandPool {
public:
    void append_state_hash(detail::StateHash&) const;
    struct Slot {
        Island island;
        u32 generation = 1;
        bool active = false;
        tick_t sleep_time_tu = 0;
        u32 toi_truncations = 0;
        tick_t frozen_until = 0;
        tick_t last_frozen_until = 0;
        u64 toi_epoch = 0;
        u8 period_bucket = 255;
        std::size_t period_position = 0;
    };

    IslandId create();
    // A nonempty island requires the paired BodyPool for destruction.
    bool destroy(IslandId id);
    bool destroy(IslandId id, BodyPool& bodies);

    const Island* get(IslandId id) const;
    void get_ids(std::vector<IslandId>& out) const;
    bool set_time(IslandId id, tick_t time);
    tick_t sleep_time(IslandId id) const;
    bool set_sleep_time(IslandId id, tick_t time);
    u32 toi_truncations(IslandId id) const;
    bool set_toi_truncations(IslandId id, u32 count);
    bool increment_toi_truncations(IslandId id);
    tick_t frozen_until(IslandId id) const;
    tick_t last_frozen_until(IslandId id) const;
    bool set_last_frozen_until(IslandId id, tick_t time);
    bool set_frozen_until(IslandId id, tick_t time);
    void reset_toi_truncations();
    bool set_schedule_location(IslandId id, u8 bucket, std::size_t position);
    bool get_schedule_location(IslandId id, u8& bucket, std::size_t& position) const;
    bool is_sleeping(IslandId id) const;
    bool has_body(IslandId id, BodyHandle body) const;
    // Reassigns a body from its previous island when necessary.
    bool add_body(IslandId id, BodyHandle body, BodyPool& bodies);
    bool remove_body(IslandId id, BodyHandle body, BodyPool& bodies);

private:
    std::vector<Slot> slots_;
    std::vector<u32> free_indices_;
    u64 toi_epoch_ = 1;
};

} // namespace sqp
