// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

// A dense pool used to store data indexed by body handles.

#pragma once
#include <array>
#include <limits>
#include <span>
#include <vector>

#include "shape/slot_table.h"
#include "sqphysics2d/detail/data_storage.h"

namespace sqp::detail {

inline constexpr u32 INVALID_BODY_DENSE_INDEX = std::numeric_limits<u32>::max();

struct BodySparseSlot {
    u32 dense_index = INVALID_BODY_DENSE_INDEX;
    u32 generation = 1;
    IslandId island;
};

// Indexing
struct BodySlots {
    std::vector<u32> dense_to_sparse;
    std::vector<BodySparseSlot> sparse;
    std::vector<u32> free_list;

    void reserve(u32 capacity);
    BodyHandle insert(BodyConfigs::Kind kind);
    u32 find(BodyHandle handle, BodyConfigs::Kind kind) const {
        if (handle.kind() != kind || handle.index() >= sparse.size())
            return INVALID_BODY_DENSE_INDEX;
        const auto& slot = sparse[handle.index()];
        return slot.generation == handle.generation ? slot.dense_index : INVALID_BODY_DENSE_INDEX;
    }
    void erase(BodyHandle handle, u32 dense_index, u32 last_dense);
    void clear();
};

struct ShapeAttachmentRange {
    u32 offset = 0;
    u32 size = 0;
    u32 capacity = 0;
};

class ShapeAttachmentArena {
public:
    void reserve(u32 capacity) {
        values_.reserve(capacity);
    }
    void clear();
    void release(ShapeAttachmentRange& range);
    void push_back(ShapeAttachmentRange& range, ShapeHandle shape);
    void erase(ShapeAttachmentRange& range, u32 index);
    std::span<ShapeHandle> get(ShapeAttachmentRange range) {
        return range.size ? std::span<ShapeHandle>{values_.data() + range.offset, range.size}
                          : std::span<ShapeHandle>{};
    }
    std::span<const ShapeHandle> get(ShapeAttachmentRange range) const {
        return range.size ? std::span<const ShapeHandle>{values_.data() + range.offset, range.size}
                          : std::span<const ShapeHandle>{};
    }

private:
    static u32 size_class(u32 capacity);
    ShapeAttachmentRange allocate(u32 capacity);
    std::vector<ShapeHandle> values_;
    std::array<std::vector<u32>, 32> free_offsets_;
};

struct RigidStorage {
    BodySlots slots;

    std::vector<real> positions_x;
    std::vector<real> positions_y;
    std::vector<real> rotations;

    std::vector<real> linear_velocities_x;
    std::vector<real> linear_velocities_y;
    std::vector<real> angular_velocities;
    std::vector<real> forces_x, forces_y, torques;

    std::vector<BodyConfigs> configs;
    std::vector<u64> motion_versions;

    std::vector<real> masses;
    std::vector<real> inertias;
    std::vector<real> inv_masses;
    std::vector<real> inv_inertias;

    std::vector<real> drag_multipliers;
    std::vector<real> local_centers_x;
    std::vector<real> local_centers_y;
    std::vector<real> bound_radii;

    std::vector<u64> category_bits;
    std::vector<u64> collision_masks;

    std::vector<ShapeAttachmentRange> shapes;
    ShapeAttachmentArena shape_arena;
    std::vector<u32> joint_heads;
    std::vector<u32> joint_counts;

    void reserve(u32 capacity);
    void erase(BodyHandle handle, u32 dense_index);
    void clear();
    u32 size() const {
        return static_cast<u32>(positions_x.size());
    }
};

// Lite bodies are the same except they only have one shape and they never rotate.

struct LiteStorage {
    BodyConfigs::Kind kind = BodyConfigs::Kind::None;
    BodySlots slots;

    std::vector<real> positions_x;
    std::vector<real> positions_y;
    std::vector<real> linear_velocities_x;
    std::vector<real> linear_velocities_y;
    std::vector<real> forces_x, forces_y;

    std::vector<BodyConfigs> configs;
    std::vector<u64> motion_versions;

    std::vector<real> masses;
    std::vector<real> inv_masses;

    std::vector<real> drag_multipliers;

    std::vector<u64> category_bits;
    std::vector<u64> collision_masks;

    std::vector<real> densities;

    std::vector<real> radii;
    std::vector<real> half_extents_x;
    std::vector<real> half_extents_y;

    // One ShapePool-owned shape per lite body; geometry below is a refreshed cache.
    std::vector<ShapeHandle> shapes;
    std::vector<u32> joint_heads;
    std::vector<u32> joint_counts;

    explicit LiteStorage(BodyConfigs::Kind value) : kind(value) {}
    void reserve(u32 capacity);
    void erase(BodyHandle handle, u32 dense_index);
    void clear();
    u32 size() const {
        return static_cast<u32>(positions_x.size());
    }
};

struct BodyStorage {
    RigidStorage rigid;
    LiteStorage circles{BodyConfigs::Kind::LiteCircle};
    LiteStorage boxes{BodyConfigs::Kind::AabbBox};
    ShapeSlotTable<BodyHandle> shape_owners;
};

} // namespace sqp::detail

namespace sqp {

struct BodyPool::Impl : detail::BodyStorage {};
} // namespace sqp
