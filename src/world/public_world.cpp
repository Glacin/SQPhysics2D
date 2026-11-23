// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "sqphysics2d/world.h"

#include "world_stages.h"
#include "state_hash.h"

#include <algorithm>
#include <cstdio>
#include <limits>
#include <utility>

namespace sqp {

namespace {
struct OperationGuard {
    bool& active;
    ~OperationGuard() { active = false; }
};
} // namespace

// All public APIs implemented here, I think people won't want to use "command" system.

World::World() : core_(new CoreWorld()), state_(&core_->world_state_) {}
World::~World() = default;
World::World(World&&) noexcept = default;
World& World::operator=(World&&) noexcept = default;

WorldQuery World::query() const noexcept {
    return WorldQuery(*state_);
}
bool World::step(tick_t dt_tu) {
    if (operation_active_)
        return false;
    operation_active_ = true;
    OperationGuard guard{operation_active_};
    snapshot_safe_ = false;
    state_->last_error_message.clear();
    state_->stage_error_detail.clear();
    last_step_tick_hashes_.clear();
    if (dt_tu > std::numeric_limits<tick_t>::max() - state_->schedule.time) {
        state_->last_error_message = "World::step: simulation time overflow";
        std::fprintf(stderr, "%s\n", state_->last_error_message.c_str());
        return false;
    }
    const tick_t target = state_->schedule.time + dt_tu;
    bool succeeded = true;
    while (state_->schedule.time < target) {
        const tick_t now = state_->schedule.time;
        const tick_t remaining = config::TU_PER_TICK - now % config::TU_PER_TICK;
        const bool reaches_boundary = remaining <= target - now;
        const tick_t next = reaches_boundary ? now + remaining : target;
        succeeded = core_->step_until(next);
        if (!succeeded)
            break;
        if (reaches_boundary)
            last_step_tick_hashes_.push_back({next / config::TU_PER_TICK,
                                              detail::hash_world_state(*state_)});
    }
    if (succeeded && dt_tu == 0)
        succeeded = core_->step_until(target);
    if (!succeeded)
        last_step_tick_hashes_.clear();
    snapshot_safe_ = succeeded;
    return succeeded;
}
bool World::capture_snapshot(WorldSnapshot& out) const {
    if (operation_active_ || !snapshot_safe_)
        return false;
    try {
        auto copy = std::make_unique<WorldState>(*state_);
        out.state_.swap(copy);
        return true;
    } catch (...) {
        return false;
    }
}
bool World::restore_snapshot(const WorldSnapshot& snapshot) {
    if (operation_active_ || !snapshot.state_)
        return false;
    try {
        auto replacement = std::unique_ptr<CoreWorld>(new CoreWorld(*snapshot.state_));
        core_.swap(replacement);
        state_ = &core_->world_state_;
        snapshot_safe_ = true;
        last_step_tick_hashes_.clear();
        return true;
    } catch (...) {
        return false;
    }
}
std::string_view World::last_error_message() const noexcept {
    return state_->last_error_message;
}
tick_t World::time_tu() const noexcept {
    return state_->schedule.time;
}
u64 World::state_hash() const {
    return detail::hash_world_state(*state_);
}
std::span<const TickHash> World::last_step_tick_hashes() const noexcept {
    return last_step_tick_hashes_;
}
const BodyPool& WorldQuery::bodies() const noexcept {
    return state_->body_pool;
}
const ShapePool& WorldQuery::shapes() const noexcept {
    return state_->shape_pool;
}
const PhysicsSettings& WorldQuery::settings() const noexcept {
    return state_->settings;
}
WorldState& World::state() noexcept {
    return *state_;
}
const WorldState& World::state() const noexcept {
    return *state_;
}

CommandId World::submit_command(tick_t time, PhysicsCommand command) {
    auto& data = state_->commands;
    if (data.next_id == 0)
        return 0;
    const CommandId id = data.next_id++;
    data.heap.push_back({time, id, std::move(command)});
    std::push_heap(data.heap.begin(), data.heap.end(), WorldState::CommandData::Later{});
    return id;
}

bool World::take_command_result(CommandId id, CommandResult& out) {
    if (id == 0)
        return false;
    auto& results = state_->commands.results;
    const auto found = std::find_if(results.begin(), results.end(),
                                    [id](const auto& item) { return item.id == id; });
    if (found == results.end())
        return false;
    out = found->result;
    results.erase(found);
    return true;
}

bool World::process_selected_commands() {
    if (operation_active_)
        return false;
    operation_active_ = true;
    OperationGuard guard{operation_active_};
    snapshot_safe_ = false;
    state_->last_error_message.clear();
    state_->stage_error_detail.clear();
    const bool succeeded = core_->process_commands();
    snapshot_safe_ = succeeded;
    return succeeded;
}
std::size_t World::queued_command_count() const noexcept {
    return state_->commands.heap.size();
}

CommandId World::create_rigid(tick_t time, const RigidBodyDef& def) {
    return submit_command(time, CreateBody{std::make_unique<RigidBodyDef>(def), {}});
}
CommandId World::create_lite_circle(tick_t time, const LiteCircleBodyDef& def, ShapeHandle shape) {
    return submit_command(time, CreateBody{std::make_unique<LiteCircleBodyDef>(def), shape});
}
CommandId World::create_aabb_box(tick_t time, const AabbBoxBodyDef& def, ShapeHandle shape) {
    return submit_command(time, CreateBody{std::make_unique<AabbBoxBodyDef>(def), shape});
}
CommandId World::destroy_body(tick_t time, BodyHandle body) {
    return submit_command(time, RemoveBody{body});
}
CommandId World::create_circle_shape(tick_t time, const CircleShapeDef& def, BodyHandle body) {
    return submit_command(time, CreateShape{def, body});
}
CommandId World::create_polygon_shape(tick_t time, const PolygonShapeDef& def, BodyHandle body) {
    return submit_command(time, CreateShape{def, body});
}
CommandId World::destroy_shape(tick_t time, ShapeHandle shape) {
    return submit_command(time, RemoveShape{shape});
}
CommandId World::set_body_velocity(tick_t time, BodyHandle body, Vec2 linear, real angular) {
    return submit_command(time, SetVelocity{body, linear, angular});
}
CommandId World::set_body_pose(tick_t time, BodyHandle body, Vec2 position, real rotation) {
    return submit_command(time, SetPose{body, position, rotation});
}
CommandId World::set_body_event_flags(tick_t time, BodyHandle body, BodyConfigs::EventFlags flags) {
    return submit_command(time, SetBodyConfig{body, flags});
}
CommandId World::set_body_override_flags(tick_t time, BodyHandle body,
                                         BodyConfigs::OverrideFlags flags) {
    return submit_command(time, SetBodyConfig{body, flags});
}
CommandId World::set_body_ccd_override(tick_t time, BodyHandle body,
                                       BodyConfigs::CcdOverride mode) {
    return submit_command(time, SetBodyConfig{body, mode});
}
CommandId World::set_body_collision_filter(tick_t time, BodyHandle body, CollisionFilter filter) {
    return submit_command(time, SetCollisionFilter{body, filter});
}
CommandId World::set_constant_force(tick_t time, BodyHandle body, Vec2 force, real torque) {
    return submit_command(time, SetConstantForce{body, force, torque});
}
CommandId World::clear_forces(tick_t time, BodyHandle body) {
    return submit_command(time, ClearForces{body});
}
CommandId World::apply_impulse(tick_t time, BodyHandle body, Vec2 impulse, real angular_impulse) {
    return submit_command(time, ApplyImpulse{body, impulse, angular_impulse});
}
CommandId World::attach_shape(tick_t time, BodyHandle body, ShapeHandle shape) {
    return submit_command(time, AttachShape{body, shape});
}
CommandId World::detach_shape(tick_t time, BodyHandle body, ShapeHandle shape) {
    return submit_command(time, DetachShape{body, shape});
}
CommandId World::set_settings(tick_t time, PhysicsSettings settings) {
    return submit_command(time, SetSettings{settings});
}
CommandId World::set_pre_solve_callback(tick_t time, EventCallback callback, void* user_data) {
    return submit_command(time, SetCallback{PhysicsEventTypes::PreSolve, callback, user_data});
}
CommandId World::set_after_solve_callback(tick_t time, EventCallback callback, void* user_data) {
    return submit_command(time, SetCallback{PhysicsEventTypes::PostSolve, callback, user_data});
}

} // namespace sqp
