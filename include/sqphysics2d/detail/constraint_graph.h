// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT
#pragma once

#include "sqphysics2d/detail/contact.h"

#include <algorithm>
#include <array>
#include <optional>
#include <tuple>
#include <vector>

namespace sqp {

// Persistent contact connectivity. Generation-checked handles prevent recycled slots
// from reviving old edges. This graph currently tracks contacts only.
class ConstraintGraph {
public:
    struct ContactEdge {
        BodyHandle body_a, body_b;
        ShapeHandle shape_a, shape_b;
        bool operator==(const ContactEdge&) const = default;
    };

    const std::vector<ContactEdge>& contacts() const noexcept {
        return contacts_;
    }
    u64 version() const noexcept { return version_; }
    const std::vector<ContactEdge>& changes() const noexcept { return changes_; }
    void clear_changes() { changes_.clear(); }
    std::vector<ContactEdge> contacts_for(BodyHandle body) const {
        std::vector<ContactEdge> result;
        if (const auto* slot = adjacency_for(body)) {
            result.reserve(slot->indices.size());
            for (std::size_t index : slot->indices)
                result.push_back(contacts_[index]);
        }
        return result;
    }

    void set_contact(const Pair& pair, bool connected) {
        ContactEdge key{pair.body_a, pair.body_b, pair.shape_a, pair.shape_b};
        const auto identity = [](const ContactEdge& edge) {
            return std::tuple{
                static_cast<u8>(edge.body_a.kind()), edge.body_a.index(), edge.body_a.generation,
                static_cast<u8>(edge.shape_a.kind),  edge.shape_a.index,  edge.shape_a.generation};
        };
        if (identity(key) > std::tuple{static_cast<u8>(key.body_b.kind()), key.body_b.index(),
                                       key.body_b.generation, static_cast<u8>(key.shape_b.kind),
                                       key.shape_b.index, key.shape_b.generation}) {
            std::swap(key.body_a, key.body_b);
            std::swap(key.shape_a, key.shape_b);
        }
        const auto index = index_of(key);
        if (connected && !index) {
            const std::size_t added = contacts_.size();
            contacts_.push_back(key);
            add_adjacency(key.body_a, added);
            if (key.body_b != key.body_a) add_adjacency(key.body_b, added);
            changes_.push_back(key);
            ++version_;
        } else if (!connected && index) {
            erase_contact(*index);
        }
    }

    template <class Predicate> void erase_if(Predicate predicate) {
        for (std::size_t index = contacts_.size(); index-- > 0;) {
            if (!predicate(contacts_[index]))
                continue;
            erase_contact(index);
        }
    }

private:
    struct AdjacencySlot {
        u32 generation = 0;
        std::vector<std::size_t> indices;
    };
    static constexpr std::size_t kind_count =
        static_cast<std::size_t>(BodyConfigs::Kind::AabbBox) + 1;

    AdjacencySlot* adjacency_for(BodyHandle body) {
        const auto kind = static_cast<std::size_t>(body.kind());
        if (kind >= adjacency_.size() || body.index() >= adjacency_[kind].size()) return nullptr;
        auto& slot = adjacency_[kind][body.index()];
        return slot.generation == body.generation ? &slot : nullptr;
    }
    const AdjacencySlot* adjacency_for(BodyHandle body) const {
        const auto kind = static_cast<std::size_t>(body.kind());
        if (kind >= adjacency_.size() || body.index() >= adjacency_[kind].size()) return nullptr;
        const auto& slot = adjacency_[kind][body.index()];
        return slot.generation == body.generation ? &slot : nullptr;
    }
    void add_adjacency(BodyHandle body, std::size_t index) {
        auto& slots = adjacency_[static_cast<std::size_t>(body.kind())];
        if (body.index() >= slots.size()) slots.resize(std::size_t(body.index()) + 1);
        auto& slot = slots[body.index()];
        if (slot.generation != body.generation) {
            slot.indices.clear();
            slot.generation = body.generation;
        }
        slot.indices.push_back(index);
    }
    std::optional<std::size_t> index_of(const ContactEdge& edge) const {
        if (const auto* slot = adjacency_for(edge.body_a))
            for (std::size_t index : slot->indices)
                if (contacts_[index] == edge) return index;
        return std::nullopt;
    }
    void forget_adjacency(BodyHandle body, std::size_t index) {
        if (auto* slot = adjacency_for(body)) std::erase(slot->indices, index);
    }
    void replace_adjacency(BodyHandle body, std::size_t from, std::size_t to) {
        if (auto* slot = adjacency_for(body))
            for (auto& index : slot->indices)
                if (index == from) {
                    index = to;
                    break;
                }
    }
    void erase_contact(std::size_t index) {
        const ContactEdge removed = contacts_[index];
        forget_adjacency(removed.body_a, index);
        if (removed.body_b != removed.body_a) forget_adjacency(removed.body_b, index);
        const std::size_t last = contacts_.size() - 1;
        if (index != last) {
            const ContactEdge moved = contacts_[last];
            contacts_[index] = moved;
            replace_adjacency(moved.body_a, last, index);
            if (moved.body_b != moved.body_a) replace_adjacency(moved.body_b, last, index);
        }
        contacts_.pop_back();
        changes_.push_back(removed);
        ++version_;
    }
    std::vector<ContactEdge> contacts_;
    std::array<std::vector<AdjacencySlot>, kind_count> adjacency_;
    std::vector<ContactEdge> changes_;
    u64 version_ = 0;
};

} // namespace sqp
