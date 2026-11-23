// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <optional>
#include <utility>
#include <vector>

#include "sqphysics2d/shape.h"

namespace sqp::detail {

template <class Value> class ShapeSlotTable {
public:
    Value* find(ShapeHandle shape) {
        auto* slot = find_slot(shape);
        return slot ? &*slot->value : nullptr;
    }
    const Value* find(ShapeHandle shape) const {
        const auto* slot = find_slot(shape);
        return slot ? &*slot->value : nullptr;
    }
    bool contains(ShapeHandle shape) const {
        return find(shape) != nullptr;
    }
    Value& at(ShapeHandle shape) {
        return *find(shape);
    }
    const Value& at(ShapeHandle shape) const {
        return *find(shape);
    }
    Value& operator[](ShapeHandle shape) {
        return try_emplace(shape).first;
    }

    template <class... Args>
    std::pair<Value&, bool> try_emplace(ShapeHandle shape, Args&&... args) {
        auto& slots = tables_[static_cast<std::size_t>(shape.kind)];
        if (slots.size() <= shape.index)
            slots.resize(static_cast<std::size_t>(shape.index) + 1);
        auto& slot = slots[shape.index];
        if (slot.generation == shape.generation && slot.value)
            return {*slot.value, false};
        slot.generation = shape.generation;
        slot.value.emplace(std::forward<Args>(args)...);
        return {*slot.value, true};
    }
    template <class... Args> std::pair<Value&, bool> emplace(ShapeHandle shape, Args&&... args) {
        return try_emplace(shape, std::forward<Args>(args)...);
    }
    void erase(ShapeHandle shape) {
        if (auto* slot = find_slot(shape))
            slot->value.reset();
    }
    void clear() {
        for (auto& table : tables_)
            for (auto& slot : table)
                slot.value.reset();
    }
    void reserve(ShapeKind kind, std::size_t size) {
        tables_[static_cast<std::size_t>(kind)].reserve(size);
    }

private:
    struct Slot {
        u32 generation = 0;
        std::optional<Value> value;
    };

    Slot* find_slot(ShapeHandle shape) {
        const auto kind = static_cast<std::size_t>(shape.kind);
        if (kind >= tables_.size() || shape.index >= tables_[kind].size())
            return nullptr;
        auto& slot = tables_[kind][shape.index];
        return slot.generation == shape.generation && slot.value ? &slot : nullptr;
    }

    const Slot* find_slot(ShapeHandle shape) const {
        const auto kind = static_cast<std::size_t>(shape.kind);
        if (kind >= tables_.size() || shape.index >= tables_[kind].size())
            return nullptr;
        const auto& slot = tables_[kind][shape.index];
        return slot.generation == shape.generation && slot.value ? &slot : nullptr;
    }

    std::array<std::vector<Slot>, 3> tables_;
};

} // namespace sqp::detail
