// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>

#include "body_storage.h"
#include "world/state_hash.h"

#include "sqphysics2d/detail/data_storage.h"

namespace sqp::detail {
namespace {
bool has_flag(BodyConfigs::OverrideFlags flags, BodyConfigs::OverrideFlags flag) {
    return (flags & flag) == flag;
}
template <class Storage> void filter_linear(Storage& s, u32 i) {
    if (has_flag(s.configs[i].override_flags, BodyConfigs::OverrideFlags::LockLinearX))
        s.linear_velocities_x[i] = real(0);
    if (has_flag(s.configs[i].override_flags, BodyConfigs::OverrideFlags::LockLinearY))
        s.linear_velocities_y[i] = real(0);
}
void update_rigid_inverse(RigidStorage& s, u32 i, bool invert_mass) {
    const bool movable =
        s.configs[i].register_type == BodyConfigs::RegisterType::Dynamic &&
        !has_flag(s.configs[i].override_flags, BodyConfigs::OverrideFlags::LockAll);
    s.inv_masses[i] = movable && invert_mass ? real(1) / s.masses[i] : real(0);
    s.inv_inertias[i] = movable && s.inertias[i] > real(0) ? real(1) / s.inertias[i] : real(0);
}
void update_lite_inverse(LiteStorage& s, u32 i) {
    s.inv_masses[i] =
        !has_flag(s.configs[i].override_flags, BodyConfigs::OverrideFlags::LockLinear) &&
                s.masses[i] > real(0)
            ? real(1) / s.masses[i]
            : real(0);
}
template <class Storage, class Callback>
bool visit_storage(Storage& storage, BodyHandle h, BodyConfigs::Kind kind, Callback&& callback) {
    const u32 i = storage.slots.find(h, kind);
    return i != INVALID_BODY_DENSE_INDEX && callback(storage, i);
}
template <class Storage, class Callback>
bool visit_body(Storage& storage, BodyHandle h, Callback&& callback) {
    switch (h.kind()) {
    case BodyConfigs::Kind::Rigid:
        return visit_storage(storage.rigid, h, h.kind(), std::forward<Callback>(callback));
    case BodyConfigs::Kind::LiteCircle:
        return visit_storage(storage.circles, h, h.kind(), std::forward<Callback>(callback));
    case BodyConfigs::Kind::AabbBox:
        return visit_storage(storage.boxes, h, h.kind(), std::forward<Callback>(callback));
    default:
        return false;
    }
}
template <class Storage>
auto select_scalar_column(Storage& s, BodyScalarColumn field)
    -> std::span<std::remove_reference_t<decltype(s.positions_x[0])>> {
    switch (field) {
    case BodyScalarColumn::PositionX:
        return s.positions_x;
    case BodyScalarColumn::PositionY:
        return s.positions_y;
    case BodyScalarColumn::VelocityX:
        return s.linear_velocities_x;
    case BodyScalarColumn::VelocityY:
        return s.linear_velocities_y;
    case BodyScalarColumn::ForceX:
        return s.forces_x;
    case BodyScalarColumn::ForceY:
        return s.forces_y;
    case BodyScalarColumn::Mass:
        return s.masses;
    case BodyScalarColumn::InverseMass:
        return s.inv_masses;
    case BodyScalarColumn::DragMultiplier:
        return s.drag_multipliers;
    default:
        break;
    }
    if constexpr (std::is_same_v<std::remove_cv_t<Storage>, RigidStorage>) {
        switch (field) {
        case BodyScalarColumn::Rotation:
            return s.rotations;
        case BodyScalarColumn::AngularVelocity:
            return s.angular_velocities;
        case BodyScalarColumn::Torque:
            return s.torques;
        case BodyScalarColumn::Inertia:
            return s.inertias;
        case BodyScalarColumn::InverseInertia:
            return s.inv_inertias;
        case BodyScalarColumn::LocalCenterX:
            return s.local_centers_x;
        case BodyScalarColumn::LocalCenterY:
            return s.local_centers_y;
        case BodyScalarColumn::BoundRadius:
            return s.bound_radii;
        default:
            break;
        }
    } else {
        switch (field) {
        case BodyScalarColumn::Density:
            return s.densities;
        case BodyScalarColumn::Radius:
            return s.radii;
        case BodyScalarColumn::HalfExtentX:
            return s.half_extents_x;
        case BodyScalarColumn::HalfExtentY:
            return s.half_extents_y;
        default:
            break;
        }
    }
    return {};
}
template <class... Vectors> void reserve_all(u32 capacity, Vectors&... values) {
    (values.reserve(capacity), ...);
}
template <class First, class... Rest> void swap_pop_all(u32 dense, First& first, Rest&... rest) {
    if (dense + 1 != first.size()) {
        first[dense] = std::move(first.back());
        ((rest[dense] = std::move(rest.back())), ...);
    }
    first.pop_back();
    (rest.pop_back(), ...);
}
template <class... Vectors> void clear_all(Vectors&... values) {
    (values.clear(), ...);
}
} // namespace

void BodySlots::reserve(u32 capacity) {
    dense_to_sparse.reserve(capacity);
    sparse.reserve(capacity);
    free_list.reserve(capacity);
}

BodyHandle BodySlots::insert(BodyConfigs::Kind kind) {
    constexpr u32 max_body_slots = u32{1} << 28;
    if (dense_to_sparse.size() >= max_body_slots)
        return {};
    u32 sparse_index;
    if (free_list.empty()) {
        if (sparse.size() >= max_body_slots)
            return {};
        sparse_index = static_cast<u32>(sparse.size());
        sparse.push_back({});
    } else {
        sparse_index = free_list.back();
        free_list.pop_back();
    }
    const u32 dense_index = static_cast<u32>(dense_to_sparse.size());
    dense_to_sparse.push_back(sparse_index);
    sparse[sparse_index].dense_index = dense_index;
    sparse[sparse_index].island = {};
    return {sparse_index, sparse[sparse_index].generation, kind};
}

void BodySlots::erase(BodyHandle handle, u32 dense_index, u32 last_dense) {
    assert(handle.index() < sparse.size() && sparse[handle.index()].dense_index == dense_index);
    if (dense_index != last_dense) {
        const u32 moved_sparse = dense_to_sparse[last_dense];
        dense_to_sparse[dense_index] = moved_sparse;
        sparse[moved_sparse].dense_index = dense_index;
    }
    dense_to_sparse.pop_back();
    BodySparseSlot& slot = sparse[handle.index()];
    slot.dense_index = INVALID_BODY_DENSE_INDEX;
    slot.island = {};
    ++slot.generation;
    if (slot.generation == 0)
        slot.generation = 1;
    free_list.push_back(handle.index());
}

void BodySlots::clear() {
    for (u32 sparse_index : dense_to_sparse) {
        BodySparseSlot& slot = sparse[sparse_index];
        slot.dense_index = INVALID_BODY_DENSE_INDEX;
        slot.island = {};
        ++slot.generation;
        if (slot.generation == 0)
            slot.generation = 1;
    }
    dense_to_sparse.clear();
    free_list.clear();
    for (u32 i = 0; i < static_cast<u32>(sparse.size()); ++i)
        free_list.push_back(i);
}

u32 ShapeAttachmentArena::size_class(u32 capacity) {
    u32 result = 0;
    while ((u32(1) << result) < capacity && result + 1 < 32)
        ++result;
    return result;
}

ShapeAttachmentRange ShapeAttachmentArena::allocate(u32 capacity) {
    const u32 index = size_class(std::max<u32>(1, capacity));
    const u32 actual = u32(1) << index;
    auto& free = free_offsets_[index];
    if (!free.empty()) {
        const u32 offset = free.back();
        free.pop_back();
        return {offset, 0, actual};
    }
    const u32 offset = static_cast<u32>(values_.size());
    values_.resize(values_.size() + actual);
    return {offset, 0, actual};
}

void ShapeAttachmentArena::release(ShapeAttachmentRange& range) {
    if (range.capacity)
        free_offsets_[size_class(range.capacity)].push_back(range.offset);
    range = {};
}

void ShapeAttachmentArena::push_back(ShapeAttachmentRange& range, ShapeHandle shape) {
    if (range.size == range.capacity) {
        auto replacement = allocate(range.capacity ? range.capacity * 2 : 1);
        std::copy_n(values_.begin() + range.offset, range.size,
                    values_.begin() + replacement.offset);
        replacement.size = range.size;
        release(range);
        range = replacement;
    }
    values_[range.offset + range.size++] = shape;
}

void ShapeAttachmentArena::erase(ShapeAttachmentRange& range, u32 index) {
    assert(index < range.size);
    std::move(values_.begin() + range.offset + index + 1,
              values_.begin() + range.offset + range.size, values_.begin() + range.offset + index);
    --range.size;
}

void ShapeAttachmentArena::clear() {
    values_.clear();
    for (auto& free : free_offsets_)
        free.clear();
}

void RigidStorage::reserve(u32 capacity) {
    slots.reserve(capacity);
    reserve_all(capacity, positions_x, positions_y, rotations, linear_velocities_x,
                linear_velocities_y, angular_velocities, forces_x, forces_y, torques, configs,
                motion_versions, masses, inertias,
                inv_masses, inv_inertias, drag_multipliers, local_centers_x, local_centers_y,
                bound_radii, category_bits, collision_masks, shapes, joint_heads, joint_counts);
    shape_arena.reserve(capacity);
}

void RigidStorage::erase(BodyHandle handle, u32 dense) {
    const u32 last = size() - 1;
    shape_arena.release(shapes[dense]);
    swap_pop_all(dense, positions_x, positions_y, rotations, linear_velocities_x,
                 linear_velocities_y, angular_velocities, forces_x, forces_y, torques, configs,
                 motion_versions, masses,
                 inertias, inv_masses, inv_inertias, drag_multipliers, local_centers_x,
                 local_centers_y, bound_radii, category_bits, collision_masks, shapes, joint_heads,
                 joint_counts);
    slots.erase(handle, dense, last);
}

void RigidStorage::clear() {
    shape_arena.clear();
    clear_all(positions_x, positions_y, rotations, linear_velocities_x, linear_velocities_y,
              angular_velocities, forces_x, forces_y, torques, configs, motion_versions, masses,
              inertias, inv_masses,
              inv_inertias, drag_multipliers, local_centers_x, local_centers_y, bound_radii,
              category_bits, collision_masks, shapes, joint_heads, joint_counts);
    slots.clear();
}

void LiteStorage::reserve(u32 capacity) {
    slots.reserve(capacity);
    reserve_all(capacity, positions_x, positions_y, linear_velocities_x, linear_velocities_y,
                forces_x, forces_y, configs, motion_versions, masses, inv_masses,
                drag_multipliers, category_bits,
                collision_masks, densities, radii, half_extents_x, half_extents_y, shapes,
                joint_heads, joint_counts);
}

void LiteStorage::erase(BodyHandle handle, u32 dense) {
    const u32 last = size() - 1;
    swap_pop_all(dense, positions_x, positions_y, linear_velocities_x, linear_velocities_y, configs,
                 forces_x, forces_y, motion_versions, masses, inv_masses, drag_multipliers,
                 category_bits,
                 collision_masks, densities, radii, half_extents_x, half_extents_y, shapes,
                 joint_heads, joint_counts);
    slots.erase(handle, dense, last);
}

void LiteStorage::clear() {
    clear_all(positions_x, positions_y, linear_velocities_x, linear_velocities_y, configs,
              forces_x, forces_y, motion_versions, masses, inv_masses, drag_multipliers,
              category_bits, collision_masks,
              densities, radii, half_extents_x, half_extents_y, shapes, joint_heads, joint_counts);
    slots.clear();
}

} // namespace sqp::detail

namespace sqp {

using namespace detail;
namespace {
constexpr auto LIVE_FLAGS = BodyConfigs::RuntimeFlags::IsAlive;

// Temporary validation result, never exposed as a BodyPool runtime value.
struct LiteGeometry {
    real radius = real(0);
    real half_x = real(0);
    real half_y = real(0);
    real density = real(0);
    real mass = real(0);
};

bool read_lite_geometry(BodyConfigs::Kind kind, ShapeHandle shape, const ShapePool& shapes,
                        LiteGeometry& out) {
    LiteGeometry candidate;
    real area = real(0);
    if (kind == BodyConfigs::Kind::LiteCircle) {
        Circle circle;
        if (!shapes.get_circle(shape, circle) || circle.center != Vec2{} ||
            !std::isfinite(circle.radius) || circle.radius <= real(0))
            return false;
        candidate.radius = circle.radius;
        area = math::PI * circle.radius * circle.radius;
    } else if (kind == BodyConfigs::Kind::AabbBox) {
        PolygonVertices vertices;
        if (!shapes.get_polygon(shape, vertices) || vertices.size() != 4)
            return false;
        real min_x = vertices[0].x, max_x = min_x;
        real min_y = vertices[0].y, max_y = min_y;
        for (const auto& v : vertices) {
            if (!std::isfinite(v.x) || !std::isfinite(v.y))
                return false;
            min_x = std::min(min_x, v.x);
            max_x = std::max(max_x, v.x);
            min_y = std::min(min_y, v.y);
            max_y = std::max(max_y, v.y);
        }
        if (max_x <= real(0) || max_y <= real(0) || min_x != -max_x || min_y != -max_y)
            return false;
        u32 corners = 0;
        for (const auto& v : vertices) {
            if ((v.x != min_x && v.x != max_x) || (v.y != min_y && v.y != max_y))
                return false;
            const u32 corner = (v.x == max_x ? 1u : 0u) | (v.y == max_y ? 2u : 0u);
            corners |= 1u << corner;
        }
        if (corners != 15u)
            return false;
        candidate.half_x = max_x;
        candidate.half_y = max_y;
        area = real(4) * max_x * max_y;
    } else {
        return false;
    }
    ShapeProperties properties;
    if (!shapes.get_properties(shape, properties))
        return false;
    candidate.density = properties.density;
    candidate.mass = area * candidate.density;
    out = candidate;
    return true;
}

void apply_geometry(LiteStorage& s, u32 i, const LiteGeometry& geometry) {
    s.radii[i] = geometry.radius;
    s.half_extents_x[i] = geometry.half_x;
    s.half_extents_y[i] = geometry.half_y;
    s.densities[i] = geometry.density;
    s.masses[i] = geometry.mass;
    update_lite_inverse(s, i);
}

template <class Def>
BodyHandle add_lite(BodyStorage& pool, LiteStorage& s, const Def& def, ShapeHandle shape,
                    const ShapePool& shapes) {
    LiteGeometry geometry;
    if (pool.shape_owners.contains(shape) || !read_lite_geometry(s.kind, shape, shapes, geometry))
        return {};
    const BodyHandle handle = s.slots.insert(s.kind);
    if (!sqp::is_valid(handle))
        return {};
    s.positions_x.push_back(def.position.x);
    s.positions_y.push_back(def.position.y);
    s.linear_velocities_x.push_back(def.linear_velocity.x);
    s.linear_velocities_y.push_back(def.linear_velocity.y);
    s.forces_x.push_back(real(0));
    s.forces_y.push_back(real(0));
    auto configs = def.configs;
    configs.kind = s.kind;
    configs.register_type = BodyConfigs::RegisterType::Dynamic;
    configs.runtime_flags = LIVE_FLAGS;
    s.configs.push_back(configs);
    s.motion_versions.push_back(0);
    s.masses.push_back(geometry.mass);
    s.inv_masses.push_back(real(0));
    s.drag_multipliers.push_back(std::max(real(0), def.drag_multiplier));
    s.category_bits.push_back(def.collision_filter.category_bits);
    s.collision_masks.push_back(def.collision_filter.collision_mask);
    s.densities.push_back(geometry.density);
    s.radii.push_back(geometry.radius);
    s.half_extents_x.push_back(geometry.half_x);
    s.half_extents_y.push_back(geometry.half_y);
    s.shapes.push_back(shape);
    s.joint_heads.push_back(INVALID_JOINT_EDGE);
    s.joint_counts.push_back(0);
    const u32 i = s.size() - 1;
    filter_linear(s, i);
    update_lite_inverse(s, i);
    pool.shape_owners.emplace(shape, handle);
    return handle;
}

void release_lite_shape(BodyStorage& pool, LiteStorage& s, u32 i) {
    pool.shape_owners.erase(s.shapes[i]);
    s.shapes[i] = {};
    apply_geometry(s, i, {});
}
} // namespace

BodyPool::BodyPool() : impl_(std::make_unique<Impl>()) {}
BodyPool::~BodyPool() = default;
BodyPool::BodyPool(const BodyPool& other) : impl_(std::make_unique<Impl>(*other.impl_)) {}
BodyPool& BodyPool::operator=(const BodyPool& other) {
    if (this != &other)
        impl_ = std::make_unique<Impl>(*other.impl_);
    return *this;
}
BodyPool::BodyPool(BodyPool&&) noexcept = default;
BodyPool& BodyPool::operator=(BodyPool&&) noexcept = default;

bool BodyPool::clear() {
    const auto has_assigned = [](const auto& storage) {
        for (u32 index : storage.slots.dense_to_sparse)
            if (storage.slots.sparse[index].island.valid())
                return true;
        return false;
    };
    if (has_assigned(impl_->rigid) || has_assigned(impl_->circles) || has_assigned(impl_->boxes))
        return false;
    const auto has_joints = [](const auto& storage) {
        return std::any_of(storage.joint_counts.begin(), storage.joint_counts.end(),
                           [](u32 count) { return count != 0; });
    };
    if (has_joints(impl_->rigid) || has_joints(impl_->circles) || has_joints(impl_->boxes))
        return false;
    impl_->rigid.clear();
    impl_->circles.clear();
    impl_->boxes.clear();
    impl_->shape_owners.clear();
    return true;
}

bool BodyPool::clear(IslandPool& islands) {
    std::vector<BodyHandle> handles;
    get_handles(handles);
    for (BodyHandle body : handles) {
        const auto id = island_of(body);
        if (id && !islands.has_body(*id, body))
            return false;
    }
    const auto has_joints = [](const auto& storage) {
        return std::any_of(storage.joint_counts.begin(), storage.joint_counts.end(),
                           [](u32 count) { return count != 0; });
    };
    if (has_joints(impl_->rigid) || has_joints(impl_->circles) || has_joints(impl_->boxes))
        return false;
    for (BodyHandle body : handles) {
        const auto id = island_of(body);
        if (id && !islands.remove_body(*id, body, *this))
            return false;
    }
    return clear();
}

bool BodyPool::clear(JointPool& joints) {
    std::vector<BodyHandle> handles;
    get_handles(handles);
    for (BodyHandle body : handles) {
        if (island_of(body))
            return false;
        u32 count = 0;
        for (auto joint = joints.first_attached(*this, body); joint.valid();
             joint = joints.next_attached(body, joint))
            ++count;
        const auto stored = visit_body(std::as_const(*impl_), body, [&](const auto& s, u32 i) {
            return count == s.joint_counts[i];
        });
        if (!stored)
            return false;
    }
    for (BodyHandle body : handles)
        if (!joints.destroy_attached(*this, body))
            return false;
    return clear();
}

bool BodyPool::clear(IslandPool& islands, JointPool& joints) {
    std::vector<BodyHandle> handles;
    get_handles(handles);
    for (BodyHandle body : handles) {
        const auto id = island_of(body);
        if (id && !islands.has_body(*id, body))
            return false;
        u32 count = 0;
        for (auto joint = joints.first_attached(*this, body); joint.valid();
             joint = joints.next_attached(body, joint))
            ++count;
        if (!visit_body(std::as_const(*impl_), body,
                        [&](const auto& s, u32 i) { return count == s.joint_counts[i]; }))
            return false;
    }
    for (BodyHandle body : handles)
        if (!joints.destroy_attached(*this, body))
            return false;
    return clear(islands);
}

void BodyPool::reserve(u32 rigid_capacity, u32 circle_capacity, u32 box_capacity) {
    impl_->rigid.reserve(rigid_capacity);
    impl_->circles.reserve(circle_capacity);
    impl_->boxes.reserve(box_capacity);
    impl_->shape_owners.reserve(ShapeKind::Circle, rigid_capacity + circle_capacity);
    impl_->shape_owners.reserve(ShapeKind::Polygon, rigid_capacity + box_capacity);
}

BodyHandle BodyPool::add_rigid(const RigidBodyDef& def) {
    RigidStorage& s = impl_->rigid;
    const BodyHandle handle = s.slots.insert(BodyConfigs::Kind::Rigid);
    if (!sqp::is_valid(handle))
        return {};
    s.positions_x.push_back(def.position.x);
    s.positions_y.push_back(def.position.y);
    s.rotations.push_back(def.rotation);
    s.linear_velocities_x.push_back(
        has_flag(def.configs.override_flags, BodyConfigs::OverrideFlags::LockLinearX)
            ? real(0)
            : def.linear_velocity.x);
    s.linear_velocities_y.push_back(
        has_flag(def.configs.override_flags, BodyConfigs::OverrideFlags::LockLinearY)
            ? real(0)
            : def.linear_velocity.y);
    s.angular_velocities.push_back(
        has_flag(def.configs.override_flags, BodyConfigs::OverrideFlags::LockRot)
            ? real(0)
            : def.angular_velocity);
    s.forces_x.push_back(real(0));
    s.forces_y.push_back(real(0));
    s.torques.push_back(real(0));
    auto configs = def.configs;
    configs.kind = BodyConfigs::Kind::Rigid;
    configs.runtime_flags = LIVE_FLAGS;
    s.configs.push_back(configs);
    s.motion_versions.push_back(0);
    s.masses.push_back(real(0));
    s.inertias.push_back(real(0));
    s.inv_masses.push_back(real(0));
    s.inv_inertias.push_back(real(0));
    s.drag_multipliers.push_back(std::max(real(0), def.drag_multiplier));
    s.local_centers_x.push_back(real(0));
    s.local_centers_y.push_back(real(0));
    s.bound_radii.push_back(real(0));
    s.category_bits.push_back(def.collision_filter.category_bits);
    s.collision_masks.push_back(def.collision_filter.collision_mask);
    s.shapes.emplace_back();
    s.joint_heads.push_back(INVALID_JOINT_EDGE);
    s.joint_counts.push_back(0);
    return handle;
}

BodyHandle BodyPool::add_lite_circle(const LiteCircleBodyDef& def, ShapeHandle shape,
                                     const ShapePool& shapes) {
    return add_lite(*impl_, impl_->circles, def, shape, shapes);
}

BodyHandle BodyPool::add_aabb_box(const AabbBoxBodyDef& def, ShapeHandle shape,
                                  const ShapePool& shapes) {
    return add_lite(*impl_, impl_->boxes, def, shape, shapes);
}

bool BodyPool::is_valid(BodyHandle h) const {
    return visit_body(std::as_const(*impl_), h, [](const auto&, u32) { return true; });
}

std::optional<u32> BodyPool::dense_index(BodyHandle h) const {
    if (!impl_)
        return std::nullopt;
    u32 index = INVALID_BODY_DENSE_INDEX;
    switch (h.kind()) {
    case BodyConfigs::Kind::Rigid:
        index = impl_->rigid.slots.find(h, h.kind());
        break;
    case BodyConfigs::Kind::LiteCircle:
        index = impl_->circles.slots.find(h, h.kind());
        break;
    case BodyConfigs::Kind::AabbBox:
        index = impl_->boxes.slots.find(h, h.kind());
        break;
    default:
        break;
    }
    return index == INVALID_BODY_DENSE_INDEX ? std::nullopt : std::optional<u32>(index);
}

std::optional<IslandId> BodyPool::island_of(BodyHandle h) const {
    if (!impl_)
        return std::nullopt;
    std::optional<IslandId> result;
    visit_body(std::as_const(*impl_), h, [&](const auto& storage, u32) {
        const IslandId id = storage.slots.sparse[h.index()].island;
        if (id.valid())
            result = id;
        return true;
    });
    return result;
}

bool BodyPool::assign_island(BodyHandle h, IslandId id) {
    if (!impl_)
        return false;
    return visit_body(*impl_, h, [&](auto& storage, u32) {
        storage.slots.sparse[h.index()].island = id;
        return true;
    });
}

std::span<real> BodyPool::scalar_column(BodyConfigs::Kind kind, BodyScalarColumn field) {
    if (!impl_)
        return {};
    switch (kind) {
    case BodyConfigs::Kind::Rigid:
        return select_scalar_column(impl_->rigid, field);
    case BodyConfigs::Kind::LiteCircle:
        return select_scalar_column(impl_->circles, field);
    case BodyConfigs::Kind::AabbBox:
        return select_scalar_column(impl_->boxes, field);
    default:
        return {};
    }
}

std::span<const real> BodyPool::scalar_column(BodyConfigs::Kind kind,
                                              BodyScalarColumn field) const {
    if (!impl_)
        return {};
    switch (kind) {
    case BodyConfigs::Kind::Rigid:
        return select_scalar_column(std::as_const(impl_->rigid), field);
    case BodyConfigs::Kind::LiteCircle:
        return select_scalar_column(std::as_const(impl_->circles), field);
    case BodyConfigs::Kind::AabbBox:
        return select_scalar_column(std::as_const(impl_->boxes), field);
    default:
        return {};
    }
}

std::span<BodyConfigs> BodyPool::configs_column(BodyConfigs::Kind kind) {
    if (!impl_)
        return {};
    switch (kind) {
    case BodyConfigs::Kind::Rigid:
        return impl_->rigid.configs;
    case BodyConfigs::Kind::LiteCircle:
        return impl_->circles.configs;
    case BodyConfigs::Kind::AabbBox:
        return impl_->boxes.configs;
    default:
        return {};
    }
}

std::span<const BodyConfigs> BodyPool::configs_column(BodyConfigs::Kind kind) const {
    if (!impl_)
        return {};
    switch (kind) {
    case BodyConfigs::Kind::Rigid:
        return impl_->rigid.configs;
    case BodyConfigs::Kind::LiteCircle:
        return impl_->circles.configs;
    case BodyConfigs::Kind::AabbBox:
        return impl_->boxes.configs;
    default:
        return {};
    }
}

std::span<u64> BodyPool::motion_version_column(BodyConfigs::Kind kind) {
    if (!impl_)
        return {};
    switch (kind) {
    case BodyConfigs::Kind::Rigid:
        return impl_->rigid.motion_versions;
    case BodyConfigs::Kind::LiteCircle:
        return impl_->circles.motion_versions;
    case BodyConfigs::Kind::AabbBox:
        return impl_->boxes.motion_versions;
    default:
        return {};
    }
}

std::span<const u64> BodyPool::motion_version_column(BodyConfigs::Kind kind) const {
    if (!impl_)
        return {};
    switch (kind) {
    case BodyConfigs::Kind::Rigid:
        return impl_->rigid.motion_versions;
    case BodyConfigs::Kind::LiteCircle:
        return impl_->circles.motion_versions;
    case BodyConfigs::Kind::AabbBox:
        return impl_->boxes.motion_versions;
    default:
        return {};
    }
}

BodyConfigs::Kind BodyPool::kind(BodyHandle h) const {
    return is_valid(h) ? h.kind() : BodyConfigs::Kind::None;
}
u32 BodyPool::rigid_count() const {
    return impl_->rigid.size();
}
u32 BodyPool::lite_circle_count() const {
    return impl_->circles.size();
}
u32 BodyPool::aabb_box_count() const {
    return impl_->boxes.size();
}
u32 BodyPool::body_count() const {
    return rigid_count() + lite_circle_count() + aabb_box_count();
}

bool BodyPool::remove(BodyHandle h) {
    if (island_of(h))
        return false;
    return visit_body(*impl_, h, [&](auto& s, u32 i) {
        if (s.joint_counts[i] != 0)
            return false;
        if constexpr (std::is_same_v<std::remove_cvref_t<decltype(s)>, RigidStorage>) {
            for (ShapeHandle shape : s.shape_arena.get(s.shapes[i]))
                impl_->shape_owners.erase(shape);
        } else {
            impl_->shape_owners.erase(s.shapes[i]);
        }
        s.erase(h, i);
        return true;
    });
}

bool BodyPool::remove(BodyHandle h, IslandPool& islands) {
    if (!is_valid(h))
        return false;
    if (visit_body(std::as_const(*impl_), h,
                   [](const auto& s, u32 i) { return s.joint_counts[i] != 0; }))
        return false;
    if (const auto id = island_of(h)) {
        if (!islands.has_body(*id, h) || !islands.remove_body(*id, h, *this))
            return false;
    }
    return remove(h);
}

bool BodyPool::remove(BodyHandle h, JointPool& joints) {
    if (!is_valid(h) || island_of(h))
        return false;
    u32 count = 0;
    for (auto joint = joints.first_attached(*this, h); joint.valid();
         joint = joints.next_attached(h, joint))
        ++count;
    if (!visit_body(std::as_const(*impl_), h,
                    [&](const auto& s, u32 i) { return count == s.joint_counts[i]; }))
        return false;
    return joints.destroy_attached(*this, h) && remove(h);
}

bool BodyPool::remove(BodyHandle h, IslandPool& islands, JointPool& joints) {
    if (!is_valid(h))
        return false;
    if (const auto id = island_of(h); id && !islands.has_body(*id, h))
        return false;
    u32 count = 0;
    for (auto joint = joints.first_attached(*this, h); joint.valid();
         joint = joints.next_attached(h, joint))
        ++count;
    if (!visit_body(std::as_const(*impl_), h,
                    [&](const auto& s, u32 i) { return count == s.joint_counts[i]; }))
        return false;
    return joints.destroy_attached(*this, h) && remove(h, islands);
}

void BodyPool::get_handles(std::vector<BodyHandle>& out) const {
    out.clear();
    out.reserve(body_count());
    const auto append = [&](const auto& storage, BodyConfigs::Kind body_kind) {
        for (u32 sparse : storage.slots.dense_to_sparse)
            out.push_back({sparse, storage.slots.sparse[sparse].generation, body_kind});
    };
    append(impl_->rigid, BodyConfigs::Kind::Rigid);
    append(impl_->circles, BodyConfigs::Kind::LiteCircle);
    append(impl_->boxes, BodyConfigs::Kind::AabbBox);
}

bool BodyPool::get_state(BodyHandle body, BodyState& out) const {
    return visit_body(std::as_const(*impl_), body, [&](const auto& s, u32 i) {
        BodyState value;
        value.position = {s.positions_x[i], s.positions_y[i]};
        value.linear_velocity = {s.linear_velocities_x[i], s.linear_velocities_y[i]};
        value.flags = s.configs[i].runtime_flags;
        if constexpr (std::same_as<std::remove_cvref_t<decltype(s)>, RigidStorage>) {
            value.rotation = s.rotations[i];
            value.angular_velocity = s.angular_velocities[i];
        }
        out = value;
        return true;
    });
}

bool BodyPool::set_state(BodyHandle body, const BodyState& value) {
    return visit_body(*impl_, body, [&](auto& s, u32 i) {
        s.positions_x[i] = value.position.x;
        s.positions_y[i] = value.position.y;
        s.linear_velocities_x[i] = value.linear_velocity.x;
        s.linear_velocities_y[i] = value.linear_velocity.y;
        s.configs[i].runtime_flags = value.flags;
        filter_linear(s, i);
        if constexpr (std::same_as<std::remove_cvref_t<decltype(s)>, RigidStorage>) {
            s.rotations[i] = value.rotation;
            s.angular_velocities[i] =
                has_flag(s.configs[i].override_flags, BodyConfigs::OverrideFlags::LockRot)
                    ? real(0)
                    : value.angular_velocity;
        }
        return true;
    });
}

bool BodyPool::get_properties(BodyHandle body, BodyProperties& out) const {
    return visit_body(std::as_const(*impl_), body, [&](const auto& s, u32 i) {
        BodyProperties value;
        value.mass = s.masses[i];
        value.inv_mass = s.inv_masses[i];
        if constexpr (std::same_as<std::remove_cvref_t<decltype(s)>, RigidStorage>) {
            value.inertia = s.inertias[i];
            value.inv_inertia = s.inv_inertias[i];
            value.local_center = {s.local_centers_x[i], s.local_centers_y[i]};
            value.bound_radius = s.bound_radii[i];
        } else {
            value.density = s.densities[i];
            value.radius = s.radii[i];
            value.half_extents = {s.half_extents_x[i], s.half_extents_y[i]};
        }
        out = value;
        return true;
    });
}

bool BodyPool::get_configuration(BodyHandle body, BodyConfiguration& out) const {
    return visit_body(std::as_const(*impl_), body, [&](const auto& s, u32 i) {
        out = {s.drag_multipliers[i], s.configs[i], {s.category_bits[i], s.collision_masks[i]}};
        return true;
    });
}

bool BodyPool::get_event_flags(BodyHandle body, BodyConfigs::EventFlags& out) const {
    return visit_body(std::as_const(*impl_), body, [&](const auto& s, u32 i) {
        out = s.configs[i].event_flags;
        return true;
    });
}

bool BodyPool::get_register_type(BodyHandle body, BodyConfigs::RegisterType& out) const {
    return visit_body(std::as_const(*impl_), body, [&](const auto& storage, u32 i) {
        if constexpr (std::same_as<std::remove_cvref_t<decltype(storage)>, RigidStorage>)
            out = storage.configs[i].register_type;
        else
            out = BodyConfigs::RegisterType::Dynamic;
        return true;
    });
}

bool BodyPool::get_constant_force(BodyHandle body, Vec2& force, real& torque) const {
    return visit_body(std::as_const(*impl_), body, [&](const auto& s, u32 i) {
        force = {s.forces_x[i], s.forces_y[i]};
        if constexpr (std::same_as<std::remove_cvref_t<decltype(s)>, RigidStorage>)
            torque = s.torques[i];
        else
            torque = real(0);
        return true;
    });
}

bool BodyPool::set_pose(BodyHandle body, Vec2 position, real rotation) {
    BodyState value;
    if (!get_state(body, value))
        return false;
    value.position = position;
    value.rotation = rotation;
    return set_state(body, value);
}
bool BodyPool::set_velocity(BodyHandle body, Vec2 linear, real angular) {
    BodyState value;
    if (!get_state(body, value))
        return false;
    value.linear_velocity = linear;
    value.angular_velocity = angular;
    return set_state(body, value);
}
bool BodyPool::set_drag_multiplier(BodyHandle body, real value) {
    return visit_body(*impl_, body, [&](auto& s, u32 i) {
        value = std::max(real(0), value);
        if (s.drag_multipliers[i] == value)
            return true;
        s.drag_multipliers[i] = value;
        return true;
    });
}
bool BodyPool::set_constant_force(BodyHandle body, Vec2 force, real torque) {
    if (!std::isfinite(force.x) || !std::isfinite(force.y) || !std::isfinite(torque))
        return false;
    return visit_body(*impl_, body, [&](auto& s, u32 i) {
        if constexpr (std::same_as<std::remove_cvref_t<decltype(s)>, RigidStorage>)
            s.torques[i] = torque;
        else if (torque != real(0))
            return false;
        s.forces_x[i] = force.x;
        s.forces_y[i] = force.y;
        return true;
    });
}
bool BodyPool::clear_forces(BodyHandle body) {
    return set_constant_force(body, {}, real(0));
}
bool BodyPool::set_ccd_override(BodyHandle body, BodyConfigs::CcdOverride value) {
    return visit_body(*impl_, body, [&](auto& s, u32 i) {
        if (s.configs[i].ccd_override == value)
            return true;
        s.configs[i].ccd_override = value;
        return true;
    });
}
bool BodyPool::set_collision_filter(BodyHandle body, CollisionFilter value) {
    return visit_body(*impl_, body, [&](auto& s, u32 i) {
        if (s.category_bits[i] == value.category_bits &&
            s.collision_masks[i] == value.collision_mask)
            return true;
        s.category_bits[i] = value.category_bits;
        s.collision_masks[i] = value.collision_mask;
        return true;
    });
}
bool BodyPool::set_event_flags(BodyHandle body, BodyConfigs::EventFlags value) {
    return visit_body(*impl_, body, [&](auto& s, u32 i) {
        s.configs[i].event_flags = value;
        return true;
    });
}
bool BodyPool::set_override_flags(BodyHandle body, BodyConfigs::OverrideFlags value) {
    return visit_body(*impl_, body, [&](auto& s, u32 i) {
        if (s.configs[i].override_flags == value)
            return true;
        s.configs[i].override_flags = value;
        filter_linear(s, i);
        if constexpr (std::same_as<std::remove_cvref_t<decltype(s)>, RigidStorage>) {
            if (has_flag(value, BodyConfigs::OverrideFlags::LockRot))
                s.angular_velocities[i] = 0;
            update_rigid_inverse(s, i, s.masses[i] > 0);
        } else {
            update_lite_inverse(s, i);
        }
        return true;
    });
}

bool BodyPool::attach_shape(BodyHandle body, ShapeHandle shape, const ShapePool& shapes) {
    if (!shapes.is_valid(shape))
        return false;
    return visit_body(*impl_, body, [&](auto& s, u32 i) {
        if constexpr (std::same_as<std::remove_cvref_t<decltype(s)>, LiteStorage>) {
            LiteGeometry geometry;
            if (!read_lite_geometry(s.kind, shape, shapes, geometry))
                return false;
            if (sqp::is_valid(s.shapes[i]))
                return s.shapes[i] == shape;
            if (impl_->shape_owners.contains(shape))
                return false;
            s.shapes[i] = shape;
            apply_geometry(s, i, geometry);
            impl_->shape_owners.emplace(shape, body);
            return true;
        } else {
            const auto* owner = impl_->shape_owners.find(shape);
            if (owner)
                return *owner == body;
            s.shape_arena.push_back(s.shapes[i], shape);
            impl_->shape_owners.emplace(shape, body);
            return refresh_properties(body, shapes);
        }
    });
}

bool BodyPool::detach_shape(BodyHandle body, ShapeHandle shape, const ShapePool& shapes) {
    return visit_body(*impl_, body, [&](auto& s, u32 i) {
        if constexpr (std::same_as<std::remove_cvref_t<decltype(s)>, LiteStorage>) {
            if (!sqp::is_valid(shape) || s.shapes[i] != shape)
                return false;
            release_lite_shape(*impl_, s, i);
            return true;
        } else {
            const auto attached = s.shape_arena.get(s.shapes[i]);
            const auto it = std::find(attached.begin(), attached.end(), shape);
            if (it == attached.end())
                return false;
            s.shape_arena.erase(s.shapes[i], static_cast<u32>(it - attached.begin()));
            impl_->shape_owners.erase(shape);
            return refresh_properties(body, shapes);
        }
    });
}

bool BodyPool::refresh_properties(BodyHandle body, const ShapePool& shapes) {
    return visit_body(*impl_, body, [&](auto& s, u32 i) {
        if constexpr (std::same_as<std::remove_cvref_t<decltype(s)>, LiteStorage>) {
            if (!sqp::is_valid(s.shapes[i]))
                return true;
            LiteGeometry geometry;
            if (!read_lite_geometry(s.kind, s.shapes[i], shapes, geometry)) {
                release_lite_shape(*impl_, s, i);
                return false;
            }
            apply_geometry(s, i, geometry);
            return true;
        } else {
            u32 attached_index = 0;
            while (attached_index < s.shapes[i].size) {
                const auto attached = s.shape_arena.get(s.shapes[i]);
                const auto shape = attached[attached_index];
                if (!shapes.is_valid(shape)) {
                    impl_->shape_owners.erase(shape);
                    s.shape_arena.erase(s.shapes[i], attached_index);
                } else {
                    impl_->shape_owners[shape] = body;
                    ++attached_index;
                }
            }

            const auto attached = s.shape_arena.get(s.shapes[i]);

            real mass = real(0);
            real geometric_bound = real(0);
            Vec2 weighted_center;
            for (ShapeHandle shape : attached) {
                ShapeProperties properties;
                if (!shapes.get_properties(shape, properties))
                    continue;
                const real shape_mass = properties.area * properties.density;
                mass += shape_mass;
                const Vec2 offset = properties.center_offset;
                weighted_center += offset * shape_mass;
                geometric_bound =
                    std::max(geometric_bound, offset.length() + properties.bound_radius);
            }

            s.masses[i] = mass;
            if (mass <= real(0)) {
                s.inertias[i] = real(0);
                s.inv_masses[i] = real(0);
                s.inv_inertias[i] = real(0);
                s.local_centers_x[i] = real(0);
                s.local_centers_y[i] = real(0);
                s.bound_radii[i] = geometric_bound;
                return true;
            }

            const Vec2 center = weighted_center / mass;
            real inertia = real(0);
            real bound = real(0);
            for (ShapeHandle shape : attached) {
                ShapeProperties properties;
                if (!shapes.get_properties(shape, properties))
                    continue;
                const Vec2 offset = properties.center_offset - center;
                inertia += properties.area * properties.density *
                           (properties.inertia_per_unit_mass + offset.length_squared());
                bound = std::max(bound, offset.length() + properties.bound_radius);
            }
            s.inertias[i] = inertia;
            s.local_centers_x[i] = center.x;
            s.local_centers_y[i] = center.y;
            s.bound_radii[i] = bound;
            // The nonpositive-mass case was handled above. Preserve the refresh path's
            // existing arithmetic for all other values, including NaN.
            update_rigid_inverse(s, i, true);
            return true;
        }
    });
}
bool BodyPool::get_attached_shapes(BodyHandle body, std::vector<ShapeHandle>& out) const {
    return visit_body(std::as_const(*impl_), body, [&](const auto& s, u32 i) {
        if constexpr (std::same_as<std::remove_cvref_t<decltype(s)>, LiteStorage>) {
            out.clear();
            if (sqp::is_valid(s.shapes[i]))
                out.push_back(s.shapes[i]);
        } else {
            const auto attached = s.shape_arena.get(s.shapes[i]);
            out.assign(attached.begin(), attached.end());
        }
        return true;
    });
}

bool BodyPool::get_shape_owner(ShapeHandle shape, BodyHandle& out) const {
    const auto* found = impl_->shape_owners.find(shape);
    if (!found || !is_valid(*found))
        return false;
    out = *found;
    return true;
}

bool BodyPool::shape_offset_from_com(BodyHandle body, ShapeHandle shape, const ShapePool& shapes,
                                     real& x, real& y) const {
    return visit_body(std::as_const(*impl_), body, [&](const auto& s, u32 i) {
        if constexpr (std::same_as<std::remove_cvref_t<decltype(s)>, LiteStorage>) {
            LiteGeometry geometry;
            if (s.shapes[i] != shape || !read_lite_geometry(s.kind, shape, shapes, geometry))
                return false;
            x = real(0);
            y = real(0);
        } else {
            const auto attached = s.shape_arena.get(s.shapes[i]);
            if (std::find(attached.begin(), attached.end(), shape) == attached.end())
                return false;
            ShapeProperties properties;
            if (!shapes.get_properties(shape, properties))
                return false;
            x = properties.center_offset.x - s.local_centers_x[i];
            y = properties.center_offset.y - s.local_centers_y[i];
        }
        return true;
    });
}

void BodyPool::append_state_hash(detail::StateHash& hash) const {
    const auto slots = [&](const detail::BodySlots& value) {
        hash.sequence(value.dense_to_sparse);
        hash.add(static_cast<u64>(value.sparse.size()));
        for (const auto& slot : value.sparse) {
            hash.add(slot.dense_index); hash.add(slot.generation);
            if (slot.dense_index != detail::INVALID_BODY_DENSE_INDEX)
                hash.add(slot.island);
        }
        hash.sequence(value.free_list);
    };
    const auto common = [&](const auto& value) {
        slots(value.slots);
        hash.sequence(value.positions_x); hash.sequence(value.positions_y);
        hash.sequence(value.linear_velocities_x); hash.sequence(value.linear_velocities_y);
        hash.sequence(value.forces_x); hash.sequence(value.forces_y);
        hash.sequence(value.configs); hash.sequence(value.motion_versions);
        hash.sequence(value.masses); hash.sequence(value.inv_masses);
        hash.sequence(value.drag_multipliers);
        hash.sequence(value.category_bits); hash.sequence(value.collision_masks);
        hash.sequence(value.joint_heads); hash.sequence(value.joint_counts);
    };
    const auto& rigid = impl_->rigid;
    common(rigid);
    hash.sequence(rigid.rotations); hash.sequence(rigid.angular_velocities);
    hash.sequence(rigid.torques); hash.sequence(rigid.inertias);
    hash.sequence(rigid.inv_inertias); hash.sequence(rigid.local_centers_x);
    hash.sequence(rigid.local_centers_y); hash.sequence(rigid.bound_radii);
    hash.add(static_cast<u64>(rigid.shapes.size()));
    for (const auto& row : rigid.shapes) {
        hash.add(row.size);
        const auto shapes = rigid.shape_arena.get(row);
        hash.sequence(shapes);
    }
    const auto lite = [&](const detail::LiteStorage& value) {
        common(value);
        hash.sequence(value.densities); hash.sequence(value.radii);
        hash.sequence(value.half_extents_x); hash.sequence(value.half_extents_y);
        hash.sequence(value.shapes);
    };
    lite(impl_->circles);
    lite(impl_->boxes);
}

} // namespace sqp
