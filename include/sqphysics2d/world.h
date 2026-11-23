// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT
#pragma once

#include "sqphysics2d/body.h"
#include "sqphysics2d/detail/data_storage.h"
#include "sqphysics2d/events.h"
#include "sqphysics2d/settings.h"
#include "sqphysics2d/shape.h"

#include <cstddef>
#include <memory>
#include <span>
#include <string_view>
#include <variant>
#include <vector>

namespace sqp {

struct WorldState;
class CoreWorld;

struct TickHash {
    tick_t tick = 0;
    u64 hash = 0;
};

// Same-process state copy. Callback user_data must remain valid after restoration.
class WorldSnapshot {
public:
    WorldSnapshot();
    ~WorldSnapshot();
    WorldSnapshot(const WorldSnapshot&);
    WorldSnapshot& operator=(const WorldSnapshot&);
    WorldSnapshot(WorldSnapshot&&) noexcept;
    WorldSnapshot& operator=(WorldSnapshot&&) noexcept;

private:
    friend class World;
    std::unique_ptr<WorldState> state_;
};

// Command payloads are applied at their requested absolute simulation time.
// A command whose time has passed when it is processed expires without changing the world.
// Set the body's linear velocity and, for rigid bodies, angular velocity (rad/s).
struct SetVelocity {
    BodyHandle body;
    Vec2 linear;
    real angular = 0;
};
// Teleport a body to a world-space pose and wake it.
struct SetPose {
    BodyHandle body;
    Vec2 position;
    real rotation = 0;
};
// Remove a body and invalidate its handle.
struct RemoveBody {
    BodyHandle body;
};
// Apply an instantaneous linear impulse and angular impulse to a dynamic body.
struct ApplyImpulse {
    BodyHandle body;
    Vec2 impulse;
    real angular_impulse = 0;
};
// Replace the persistent force and torque acting on a body.
struct SetConstantForce {
    BodyHandle body;
    Vec2 force;
    real torque = 0;
};
// Remove a body's persistent force and torque.
struct ClearForces {
    BodyHandle body;
};
// Replace the body's mutual category/mask collision filter.
struct SetCollisionFilter {
    BodyHandle body;
    CollisionFilter filter;
};
// Change CCD mode, event flags, or physics override flags.
struct SetBodyConfig {
    using Value = std::variant<BodyConfigs::CcdOverride, BodyConfigs::EventFlags,
                               BodyConfigs::OverrideFlags>;
    BodyHandle body;
    Value value = BodyConfigs::CcdOverride::Auto;
};
// Attach or detach a shape; only rigid bodies support these commands.
struct AttachShape {
    BodyHandle body;
    ShapeHandle shape;
};
struct DetachShape {
    BodyHandle body;
    ShapeHandle shape;
};
// Create a body. Lite bodies require a compatible existing shape; rigid bodies start empty.
struct CreateBody {
    std::variant<std::unique_ptr<RigidBodyDef>, std::unique_ptr<LiteCircleBodyDef>,
                 std::unique_ptr<AabbBoxBodyDef>>
        definition;
    ShapeHandle shape; // Required for lite bodies; ignored for rigid bodies.
};
// Create a shape, optionally attaching it to a rigid body.
struct CreateShape {
    std::variant<CircleShapeDef, PolygonShapeDef> definition;
    BodyHandle body; // Invalid means create detached.
};
// Remove a shape, detaching it from its owner first if necessary.
struct RemoveShape {
    ShapeHandle shape;
};
// Replace world settings; the requested settings are sanitized when applied.
struct SetSettings {
    PhysicsSettings settings;
};
// Install or clear a PreSolve or PostSolve callback. A null callback clears it.
struct SetCallback {
    PhysicsEventTypes type = PhysicsEventTypes::PreSolve;
    EventCallback callback = nullptr;
    void* user_data = nullptr;
};

// Closed set of state changes accepted by World::submit_command().
using PhysicsCommand =
    std::variant<SetVelocity, SetPose, SetBodyConfig, RemoveBody, ApplyImpulse, SetConstantForce,
                 ClearForces, SetCollisionFilter, AttachShape, DetachShape, CreateBody, CreateShape,
                 RemoveShape, SetSettings, SetCallback>;

using CommandId = u64;
// Applied: the command succeeded; Failed: validation/application failed;
// Expired: its requested time passed before it was processed.
enum class CommandStatus : u8 { Applied, Failed, Expired };
// One-shot result. body/shape is populated by a successful creation command.
struct CommandResult {
    CommandStatus status = CommandStatus::Failed;
    BodyHandle body;
    ShapeHandle shape;
};

// Read-only view of the world's live pools and settings. References are not snapshots.
class WorldQuery {
public:
    // Pool references remain tied to this world and reflect later mutations.
    const BodyPool& bodies() const noexcept;
    const ShapePool& shapes() const noexcept;
    const PhysicsSettings& settings() const noexcept;

private:
    friend class World;
    explicit WorldQuery(const WorldState& state) noexcept : state_(&state) {}
    const WorldState* state_;
};

// Owns the simulation and its timestamped command queue. Move-only.
// Body and shape handles must be checked against the destination world's pools after mutations.
class World {
public:
    World();
    ~World();
    World(World&&) noexcept;
    World& operator=(World&&) noexcept;
    World(const World&) = delete;
    World& operator=(const World&) = delete;

    // Inspect the current live state without mutating it through this view.
    WorldQuery query() const noexcept;
    // Advance by dt simulation time units. Returns false if time overflows or a stage fails.
    // Queries expose each body's last solved state.
    bool step(tick_t dt_tu);
    // Capture only between successful public calls; false leaves out unchanged.
    bool capture_snapshot(WorldSnapshot& out) const;
    // Restoring invalidates previously obtained queries, pool references and spans.
    // False leaves this world unchanged. Neither method may run from a callback.
    bool restore_snapshot(const WorldSnapshot& snapshot);
    // Diagnostic for the most recent failed step or command-processing call.
    // The view remains valid until the next such call or destruction of this world.
    std::string_view last_error_message() const noexcept;
    // Current absolute simulation time in TU.
    tick_t time_tu() const noexcept;
    // Hash of the current state. Results are comparable within the same build/process.
    u64 state_hash() const;
    // One entry per complete tick boundary reached by the last successful step().
    // The span is invalidated by the next step, restore, move, or destruction.
    std::span<const TickHash> last_step_tick_hashes() const noexcept;
    // Queue a command for an absolute TU. Zero means no command ID was available.
    CommandId submit_command(tick_t time, PhysicsCommand command);
    // Consume a completed result once; false means the ID has no available result.
    bool take_command_result(CommandId id, CommandResult& out);
    // Process commands at the selected interval time. Used by the stage pipeline.
    bool process_selected_commands();
    // Number of commands still waiting in the queue.
    std::size_t queued_command_count() const noexcept;

    // Creation helpers return command IDs, not handles. Retrieve handles from CommandResult.
    // Rigid bodies start without shapes; lite bodies require the supplied compatible shape.
    CommandId create_rigid(tick_t time, const RigidBodyDef& def = {});
    CommandId create_lite_circle(tick_t time, const LiteCircleBodyDef& def, ShapeHandle shape);
    CommandId create_aabb_box(tick_t time, const AabbBoxBodyDef& def, ShapeHandle shape);
    // Destroy the named body when the command is applied.
    CommandId destroy_body(tick_t time, BodyHandle body);
    // An invalid body handle creates a detached shape.
    CommandId create_circle_shape(tick_t time, const CircleShapeDef& def, BodyHandle body = {});
    CommandId create_polygon_shape(tick_t time, const PolygonShapeDef& def, BodyHandle body = {});
    // Destroy the shape and detach it from its owner if attached.
    CommandId destroy_shape(tick_t time, ShapeHandle shape);
    // Change world-space motion or pose at the requested absolute TU.
    CommandId set_body_velocity(tick_t time, BodyHandle body, Vec2 linear, real angular = 0);
    CommandId set_body_pose(tick_t time, BodyHandle body, Vec2 position, real rotation = 0);
    // Configure event delivery, motion overrides, CCD, or collision filtering.
    CommandId set_body_event_flags(tick_t time, BodyHandle body, BodyConfigs::EventFlags flags);
    CommandId set_body_override_flags(tick_t time, BodyHandle body,
                                      BodyConfigs::OverrideFlags flags);
    CommandId set_body_ccd_override(tick_t time, BodyHandle body, BodyConfigs::CcdOverride mode);
    CommandId set_body_collision_filter(tick_t time, BodyHandle body, CollisionFilter filter);
    // Persistent force/torque lasts until replaced or cleared; impulses act immediately.
    CommandId set_constant_force(tick_t time, BodyHandle body, Vec2 force, real torque = 0);
    CommandId clear_forces(tick_t time, BodyHandle body);
    CommandId apply_impulse(tick_t time, BodyHandle body, Vec2 impulse, real angular_impulse = 0);
    // Shape attachment changes a rigid body's derived mass properties.
    CommandId attach_shape(tick_t time, BodyHandle body, ShapeHandle shape);
    CommandId detach_shape(tick_t time, BodyHandle body, ShapeHandle shape);
    // Replace settings or solver callbacks when their commands are applied.
    CommandId set_settings(tick_t time, PhysicsSettings settings);
    CommandId set_pre_solve_callback(tick_t time, EventCallback callback,
                                     void* user_data = nullptr);
    CommandId set_after_solve_callback(tick_t time, EventCallback callback,
                                       void* user_data = nullptr);

    // Direct access to stage data for advanced integration and diagnostics.
    // Mutating it bypasses World command validation.
    WorldState& state() noexcept;
    const WorldState& state() const noexcept;

private:
    std::unique_ptr<CoreWorld> core_;
    WorldState* state_ = nullptr;
    bool operation_active_ = false;
    bool snapshot_safe_ = true;
    std::vector<TickHash> last_step_tick_hashes_;
};

} // namespace sqp
