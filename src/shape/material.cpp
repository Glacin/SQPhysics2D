// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "sqphysics2d/detail/data_storage.h"

namespace sqp {

MaterialPool::MaterialPool() : materials_(1) {}

u32 MaterialPool::add(const PhysicsMaterial& value) {
    if (materials_.size() >= INVALID_MATERIAL_INDEX)
        return INVALID_MATERIAL_INDEX;
    materials_.push_back(value);
    return static_cast<u32>(materials_.size() - 1);
}

bool MaterialPool::get(u32 index, PhysicsMaterial& out) const {
    if (!contains(index))
        return false;
    out = materials_[index];
    return true;
}

bool MaterialPool::set(u32 index, const PhysicsMaterial& value) {
    if (!contains(index))
        return false;
    materials_[index] = value;
    return true;
}

bool MaterialPool::contains(u32 index) const {
    return index < materials_.size();
}

u32 MaterialPool::size() const {
    return static_cast<u32>(materials_.size());
}

void MaterialPool::reserve(u32 capacity) {
    materials_.reserve(capacity);
}

void MaterialPool::reset() {
    materials_.assign(1, PhysicsMaterial{});
}

} // namespace sqp
