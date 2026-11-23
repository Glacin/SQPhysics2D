// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "sqphysics2d/detail/data_storage.h"
#include "world/state_hash.h"

#include <algorithm>
#include <limits>

namespace sqp {
namespace {

template <class Small> auto small_end(Small& small) {
    return std::find(small.data.begin(), small.data.end(), BodyHandle{});
}

std::size_t member_count(const Island& island) {
    if (const auto* small = std::get_if<Island::SmallIsland>(&island.bodies))
        return static_cast<std::size_t>(small_end(*small) - small->data.begin());
    return std::get<std::vector<BodyHandle>>(island.bodies).size();
}

bool contains(const Island& island, BodyHandle body) {
    if (const auto* small = std::get_if<Island::SmallIsland>(&island.bodies)) {
        const auto end = small_end(*small);
        return std::find(small->data.begin(), end, body) != end;
    }
    const auto& members = std::get<std::vector<BodyHandle>>(island.bodies);
    return std::find(members.begin(), members.end(), body) != members.end();
}

void append(Island& island, BodyHandle body) {
    if (auto* small = std::get_if<Island::SmallIsland>(&island.bodies)) {
        const auto end = small_end(*small);
        if (end != small->data.end()) {
            *end = body;
            return;
        }
        std::vector<BodyHandle> members;
        members.reserve(3);
        members.insert(members.end(), small->data.begin(), small->data.end());
        members.push_back(body);
        island.bodies = std::move(members);
        return;
    }
    std::get<std::vector<BodyHandle>>(island.bodies).push_back(body);
}

void erase(Island& island, BodyHandle body) {
    if (auto* small = std::get_if<Island::SmallIsland>(&island.bodies)) {
        const auto end = small_end(*small);
        const auto found = std::find(small->data.begin(), end, body);
        std::move(found + 1, end, found);
        *(end - 1) = {};
        return;
    }
    auto& members = std::get<std::vector<BodyHandle>>(island.bodies);
    members.erase(std::find(members.begin(), members.end(), body));
    if (members.size() <= 2) {
        Island::SmallIsland small;
        std::copy(members.begin(), members.end(), small.data.begin());
        island.bodies = small;
    }
}

std::vector<BodyHandle> member_snapshot(const Island& island) {
    if (const auto* small = std::get_if<Island::SmallIsland>(&island.bodies))
        return {small->data.begin(), small_end(*small)};
    return std::get<std::vector<BodyHandle>>(island.bodies);
}

} // namespace

IslandId IslandPool::create() {
    u32 index;
    if (free_indices_.empty()) {
        if (slots_.size() >= std::numeric_limits<u32>::max())
            return {};
        index = static_cast<u32>(slots_.size());
        slots_.emplace_back();
    } else {
        index = free_indices_.back();
        free_indices_.pop_back();
    }
    Slot& slot = slots_[index];
    slot.island = {};
    slot.sleep_time_tu = 0;
    slot.toi_truncations = 0;
    slot.frozen_until = 0;
    slot.last_frozen_until = 0;
    slot.toi_epoch = toi_epoch_;
    slot.period_bucket = 255;
    slot.period_position = 0;
    slot.active = true;
    return {index, slot.generation};
}

const Island* IslandPool::get(IslandId id) const {
    if (!id.valid() || id.index >= slots_.size())
        return nullptr;
    const Slot& slot = slots_[id.index];
    return slot.active && slot.generation == id.generation ? &slot.island : nullptr;
}

void IslandPool::get_ids(std::vector<IslandId>& out) const {
    out.clear();
    for (u32 index = 0; index < slots_.size(); ++index)
        if (slots_[index].active)
            out.push_back({index, slots_[index].generation});
}

tick_t IslandPool::sleep_time(IslandId id) const {
    return get(id) ? slots_[id.index].sleep_time_tu : 0;
}

bool IslandPool::set_sleep_time(IslandId id, tick_t time) {
    if (!get(id))
        return false;
    slots_[id.index].sleep_time_tu = time;
    return true;
}

u32 IslandPool::toi_truncations(IslandId id) const {
    return get(id) && slots_[id.index].toi_epoch == toi_epoch_
               ? slots_[id.index].toi_truncations : 0;
}

bool IslandPool::set_toi_truncations(IslandId id, u32 count) {
    if (!get(id))
        return false;
    slots_[id.index].toi_truncations = count;
    slots_[id.index].toi_epoch = toi_epoch_;
    return true;
}

bool IslandPool::increment_toi_truncations(IslandId id) {
    if (!get(id))
        return false;
    if (slots_[id.index].toi_epoch != toi_epoch_) {
        slots_[id.index].toi_epoch = toi_epoch_;
        slots_[id.index].toi_truncations = 0;
    }
    if (slots_[id.index].toi_truncations == std::numeric_limits<u32>::max())
        return false;
    ++slots_[id.index].toi_truncations;
    return true;
}

tick_t IslandPool::frozen_until(IslandId id) const {
    return get(id) ? slots_[id.index].frozen_until : 0;
}

tick_t IslandPool::last_frozen_until(IslandId id) const {
    return get(id) ? slots_[id.index].last_frozen_until : 0;
}

bool IslandPool::set_last_frozen_until(IslandId id, tick_t time) {
    if (!get(id))
        return false;
    slots_[id.index].last_frozen_until = time;
    return true;
}

bool IslandPool::set_frozen_until(IslandId id, tick_t time) {
    if (!get(id))
        return false;
    slots_[id.index].frozen_until = time;
    if (time != 0)
        slots_[id.index].last_frozen_until = time;
    return true;
}

void IslandPool::reset_toi_truncations() {
    if (++toi_epoch_ == 0) {
        toi_epoch_ = 1;
        for (Slot& slot : slots_)
            slot.toi_epoch = 0;
    }
}

bool IslandPool::set_schedule_location(IslandId id, u8 bucket, std::size_t position) {
    if (!get(id))
        return false;
    Slot& slot = slots_[id.index];
    slot.period_bucket = bucket;
    slot.period_position = position;
    return true;
}

bool IslandPool::get_schedule_location(IslandId id, u8& bucket, std::size_t& position) const {
    if (!get(id))
        return false;
    const Slot& slot = slots_[id.index];
    bucket = slot.period_bucket;
    position = slot.period_position;
    return true;
}

bool IslandPool::is_sleeping(IslandId id) const {
    return get(id) && slots_[id.index].period_bucket == 0;
}

bool IslandPool::set_time(IslandId id, tick_t time) {
    if (!get(id))
        return false;
    slots_[id.index].island.time = time;
    return true;
}

bool IslandPool::has_body(IslandId id, BodyHandle body) const {
    const Island* island = get(id);
    return island && contains(*island, body);
}

bool IslandPool::add_body(IslandId id, BodyHandle body, BodyPool& bodies) {
    if (!get(id) || !bodies.is_valid(body))
        return false;
    const auto previous = bodies.island_of(body);
    if (previous == id)
        return has_body(id, body);
    if (has_body(id, body) || (previous && !has_body(*previous, body)))
        return false;

    append(slots_[id.index].island, body);
    slots_[id.index].sleep_time_tu = 0;
    if (previous)
        erase(slots_[previous->index].island, body);
    if (previous)
        slots_[previous->index].sleep_time_tu = 0;
    if (!bodies.assign_island(body, id))
        return false;
    return true;
}

bool IslandPool::remove_body(IslandId id, BodyHandle body, BodyPool& bodies) {
    if (!get(id) || !bodies.is_valid(body) || bodies.island_of(body) != id || !has_body(id, body))
        return false;
    erase(slots_[id.index].island, body);
    slots_[id.index].sleep_time_tu = 0;
    if (!bodies.assign_island(body, {}))
        return false;
    return true;
}

bool IslandPool::destroy(IslandId id) {
    if (!get(id) || member_count(slots_[id.index].island) != 0)
        return false;
    Slot& slot = slots_[id.index];
    slot.island = {};
    slot.sleep_time_tu = 0;
    slot.toi_truncations = 0;
    slot.last_frozen_until = 0;
    slot.frozen_until = 0;
    slot.period_bucket = 255;
    slot.period_position = 0;
    slot.active = false;
    if (++slot.generation == 0)
        slot.generation = 1;
    free_indices_.push_back(id.index);
    return true;
}

bool IslandPool::destroy(IslandId id, BodyPool& bodies) {
    const Island* island = get(id);
    if (!island)
        return false;
    const auto members = member_snapshot(*island);
    for (BodyHandle body : members)
        if (!bodies.is_valid(body) || bodies.island_of(body) != id)
            return false;
    for (BodyHandle body : members)
        if (!remove_body(id, body, bodies))
            return false;
    return destroy(id);
}

void IslandPool::append_state_hash(detail::StateHash& hash) const {
    hash.add(static_cast<u64>(slots_.size()));
    for (const Slot& slot : slots_) {
        hash.add(slot.generation); hash.add(slot.active);
        if (slot.active) {
            hash.add(slot.island.time);
            if (const auto* small = std::get_if<Island::SmallIsland>(&slot.island.bodies))
                hash.sequence(small->data);
            else
                hash.sequence(std::get<std::vector<BodyHandle>>(slot.island.bodies));
            hash.add(slot.sleep_time_tu); hash.add(slot.toi_truncations);
            hash.add(slot.frozen_until); hash.add(slot.last_frozen_until);
            hash.add(slot.toi_epoch); hash.add(slot.period_bucket);
            hash.add(static_cast<u64>(slot.period_position));
        }
    }
    hash.sequence(free_indices_);
    hash.add(toi_epoch_);
}

} // namespace sqp
