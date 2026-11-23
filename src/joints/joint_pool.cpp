// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "sqphysics2d/detail/data_storage.h"
#include "world/state_hash.h"

#include <cassert>
#include <utility>

#include "body/body_storage.h"

namespace sqp {
namespace {

struct BodyLinks {
    u32* head = nullptr;
    u32* count = nullptr;
};

struct ConstBodyLinks {
    const u32* head = nullptr;
    const u32* count = nullptr;
};

ConstBodyLinks links(const detail::BodyStorage& storage, BodyHandle body) {
    const auto visit = [&](const auto& s) -> ConstBodyLinks {
        const u32 dense = s.slots.find(body, body.kind());
        if (dense == detail::INVALID_BODY_DENSE_INDEX)
            return {};
        return {&s.joint_heads[dense], &s.joint_counts[dense]};
    };
    switch (body.kind()) {
    case BodyConfigs::Kind::Rigid:
        return visit(storage.rigid);
    case BodyConfigs::Kind::LiteCircle:
        return visit(storage.circles);
    case BodyConfigs::Kind::AabbBox:
        return visit(storage.boxes);
    default:
        return {};
    }
}

BodyLinks links(detail::BodyStorage& storage, BodyHandle body) {
    const auto found = links(std::as_const(storage), body);
    return {const_cast<u32*>(found.head), const_cast<u32*>(found.count)};
}

template <class T> struct DenseJoints {
    std::vector<T> values;
    std::vector<u32> owner_slots;
};

} // namespace

struct JointPool::Impl {
    struct Slot {
        JointRegistration registration;
        u32 generation = 1;
        bool active = false;
    };

    std::vector<Slot> slots;
    std::vector<u32> free_list;
    detail::BodyStorage* owner = nullptr;
    u32 live_count = 0;
    DenseJoints<PinJoint> pins;
    DenseJoints<SliderJoint> sliders;
    DenseJoints<WeldJoint> welds;
    DenseJoints<DistanceJoint> distances;
    DenseJoints<SpringJoint> springs;

    JointEdge& edge(u32 key) {
        return slots[key >> 1].registration.edges[key & 1];
    }

    JointHandle handle(u32 key) const {
        const u32 slot_index = key >> 1;
        if (key == INVALID_JOINT_EDGE || slot_index >= slots.size() || !slots[slot_index].active)
            return {};
        return {slot_index, slots[slot_index].generation};
    }

    template <class T> DenseJoints<T>& dense() {
        if constexpr (std::is_same_v<T, PinJoint>)
            return pins;
        else if constexpr (std::is_same_v<T, SliderJoint>)
            return sliders;
        else if constexpr (std::is_same_v<T, WeldJoint>)
            return welds;
        else if constexpr (std::is_same_v<T, DistanceJoint>)
            return distances;
        else
            return springs;
    }

    void link(detail::BodyStorage& bodies, BodyHandle body, u32 key) {
        const BodyLinks entry = links(bodies, body);
        assert(entry.head && entry.count);
        JointEdge& inserted = edge(key);
        inserted.prev = INVALID_JOINT_EDGE;
        inserted.next = *entry.head;
        if (*entry.head != INVALID_JOINT_EDGE)
            edge(*entry.head).prev = key;
        *entry.head = key;
        ++*entry.count;
    }

    void unlink(detail::BodyStorage& bodies, BodyHandle body, u32 key) {
        const BodyLinks entry = links(bodies, body);
        assert(entry.head && entry.count && *entry.count);
        const JointEdge old = edge(key);
        if (old.prev == INVALID_JOINT_EDGE)
            *entry.head = old.next;
        else
            edge(old.prev).next = old.next;
        if (old.next != INVALID_JOINT_EDGE)
            edge(old.next).prev = old.prev;
        edge(key) = {};
        --*entry.count;
    }

    template <class T>
    JointHandle insert(detail::BodyStorage& bodies, const JointDef& def, T runtime,
                       JointType type) {
        if (owner && owner != &bodies)
            return {};
        auto& typed = dense<T>();
        if (free_list.empty() && slots.size() > max_joint_slot)
            return {};
        const auto next_size = slots.size() + (free_list.empty() ? 1 : 0);
        slots.reserve(next_size);
        free_list.reserve(next_size);
        typed.values.reserve(typed.values.size() + 1);
        typed.owner_slots.reserve(typed.owner_slots.size() + 1);

        u32 slot_index;
        if (free_list.empty()) {
            slot_index = static_cast<u32>(slots.size());
            slots.emplace_back();
        } else {
            slot_index = free_list.back();
            free_list.pop_back();
        }
        Slot& slot = slots[slot_index];
        const u32 dense_index = static_cast<u32>(typed.values.size());
        typed.values.push_back(std::move(runtime));
        typed.owner_slots.push_back(slot_index);
        slot.registration = {def.body_a, def.body_b, {}, dense_index, type, def.collide_connected};
        slot.active = true;
        owner = &bodies;
        ++live_count;
        const JointHandle handle{slot_index, slot.generation};
        if (def.body_a != BodyHandle{})
            link(bodies, def.body_a, joint_edge_key(handle, 0));
        link(bodies, def.body_b, joint_edge_key(handle, 1));
        return handle;
    }

    void erase_type(const JointRegistration& registration) {
        u32 moved_owner = INVALID_JOINT_EDGE;
        const auto erase = [&](auto& typed) {
            // The moved owner is obtained before the last element is removed.
            const u32 last = static_cast<u32>(typed.values.size() - 1);
            if (registration.type_index != last)
                moved_owner = typed.owner_slots[last];
            if (registration.type_index != last) {
                typed.values[registration.type_index] = std::move(typed.values[last]);
                typed.owner_slots[registration.type_index] = moved_owner;
            }
            typed.values.pop_back();
            typed.owner_slots.pop_back();
        };
        switch (registration.type) {
        case JointType::Pin:
            erase(pins);
            break;
        case JointType::Slider:
            erase(sliders);
            break;
        case JointType::Weld:
            erase(welds);
            break;
        case JointType::Distance:
            erase(distances);
            break;
        case JointType::Spring:
            erase(springs);
            break;
        }
        if (moved_owner != INVALID_JOINT_EDGE)
            slots[moved_owner].registration.type_index = registration.type_index;
    }
};

JointPool::JointPool() : impl_(std::make_unique<Impl>()) {}
JointPool::~JointPool() = default;
JointPool::JointPool(JointPool&&) noexcept = default;
JointPool& JointPool::operator=(JointPool&&) noexcept = default;
void JointPool::copy_from(const JointPool& other, BodyPool& destination_bodies) {
    auto copy = std::make_unique<Impl>(*other.impl_);
    copy->owner = other.impl_->owner ? destination_bodies.impl_.get() : nullptr;
    impl_ = std::move(copy);
}

namespace {
real reference_angle(const BodyPool& bodies, const JointDef& def) {
    BodyState b;
    bodies.get_state(def.body_b, b);
    if (def.body_a == BodyHandle{})
        return b.rotation;
    BodyState a;
    bodies.get_state(def.body_a, a);
    return b.rotation - a.rotation;
}
} // namespace

JointHandle JointPool::create(BodyPool& bodies, const PinJointDef& def) {
    if (!validate_joint_def(def, bodies))
        return {};
    return impl_->insert(*bodies.impl_, def, PinJoint{def.local_anchor_a, def.local_anchor_b},
                         JointType::Pin);
}

JointHandle JointPool::create(BodyPool& bodies, const SliderJointDef& def) {
    if (!validate_joint_def(def, bodies))
        return {};
    SliderJoint runtime{def.local_anchor_a,    def.local_anchor_b,
                        def.local_axis_a,      reference_angle(bodies, def),
                        def.lower_translation, def.upper_translation,
                        def.enable_limit};
    return impl_->insert(*bodies.impl_, def, runtime, JointType::Slider);
}

JointHandle JointPool::create(BodyPool& bodies, const WeldJointDef& def) {
    if (!validate_joint_def(def, bodies))
        return {};
    WeldJoint runtime{def.local_anchor_a, def.local_anchor_b, reference_angle(bodies, def)};
    return impl_->insert(*bodies.impl_, def, runtime, JointType::Weld);
}

JointHandle JointPool::create(BodyPool& bodies, const DistanceJointDef& def) {
    if (!validate_joint_def(def, bodies))
        return {};
    DistanceJoint runtime{def.local_anchor_a, def.local_anchor_b, def.length};
    return impl_->insert(*bodies.impl_, def, runtime, JointType::Distance);
}

JointHandle JointPool::create(BodyPool& bodies, const SpringJointDef& def) {
    if (!validate_joint_def(def, bodies))
        return {};
    SpringJoint runtime{def.local_anchor_a, def.local_anchor_b, def.rest_length, def.hertz,
                        def.damping_ratio};
    return impl_->insert(*bodies.impl_, def, runtime, JointType::Spring);
}

const JointRegistration* JointPool::get(JointHandle handle) const {
    if (!handle.valid() || handle.index >= impl_->slots.size())
        return nullptr;
    const auto& slot = impl_->slots[handle.index];
    return slot.active && slot.generation == handle.generation ? &slot.registration : nullptr;
}

bool JointPool::is_valid(JointHandle handle) const {
    return get(handle) != nullptr;
}

const PinJoint* JointPool::get_pin(JointHandle h) const {
    const auto* r = get(h);
    return r && r->type == JointType::Pin ? &impl_->pins.values[r->type_index] : nullptr;
}
const SliderJoint* JointPool::get_slider(JointHandle h) const {
    const auto* r = get(h);
    return r && r->type == JointType::Slider ? &impl_->sliders.values[r->type_index] : nullptr;
}
const WeldJoint* JointPool::get_weld(JointHandle h) const {
    const auto* r = get(h);
    return r && r->type == JointType::Weld ? &impl_->welds.values[r->type_index] : nullptr;
}
const DistanceJoint* JointPool::get_distance(JointHandle h) const {
    const auto* r = get(h);
    return r && r->type == JointType::Distance ? &impl_->distances.values[r->type_index] : nullptr;
}
const SpringJoint* JointPool::get_spring(JointHandle h) const {
    const auto* r = get(h);
    return r && r->type == JointType::Spring ? &impl_->springs.values[r->type_index] : nullptr;
}

JointHandle JointPool::first_attached(const BodyPool& bodies, BodyHandle body) const {
    if (impl_->owner && impl_->owner != bodies.impl_.get())
        return {};
    const auto entry = links(*bodies.impl_, body);
    if (!entry.head)
        return {};
    const auto handle = impl_->handle(*entry.head);
    const auto* registration = get(handle);
    return registration && (registration->body_a == body || registration->body_b == body)
               ? handle
               : JointHandle{};
}

JointHandle JointPool::next_attached(BodyHandle body, JointHandle joint) const {
    const auto* registration = get(joint);
    if (!registration)
        return {};
    u32 side;
    if (registration->body_a == body)
        side = 0;
    else if (registration->body_b == body)
        side = 1;
    else
        return {};
    const auto next = impl_->handle(registration->edges[side].next);
    const auto* next_registration = get(next);
    return next_registration &&
                   (next_registration->body_a == body || next_registration->body_b == body)
               ? next
               : JointHandle{};
}

bool JointPool::destroy(BodyPool& bodies, JointHandle joint) {
    if (impl_->owner != bodies.impl_.get())
        return false;
    const auto* found = get(joint);
    if (!found)
        return false;
    const JointRegistration registration = *found;
    if (!links(*bodies.impl_, registration.body_b).head ||
        (registration.body_a != BodyHandle{} && !links(*bodies.impl_, registration.body_a).head))
        return false;
    if (registration.body_a != BodyHandle{})
        impl_->unlink(*bodies.impl_, registration.body_a, joint_edge_key(joint, 0));
    impl_->unlink(*bodies.impl_, registration.body_b, joint_edge_key(joint, 1));
    impl_->erase_type(registration);
    auto& slot = impl_->slots[joint.index];
    slot.registration = {};
    slot.active = false;
    if (++slot.generation == 0)
        slot.generation = 1;
    impl_->free_list.push_back(joint.index);
    if (--impl_->live_count == 0)
        impl_->owner = nullptr;
    return true;
}

bool JointPool::destroy_attached(BodyPool& bodies, BodyHandle body) {
    if (!bodies.is_valid(body) || (impl_->owner && impl_->owner != bodies.impl_.get()))
        return false;
    for (;;) {
        const auto entry = links(*bodies.impl_, body);
        if (*entry.head == INVALID_JOINT_EDGE)
            return true;
        const auto joint = first_attached(bodies, body);
        if (!joint.valid() || !destroy(bodies, joint))
            return false;
    }
}

void JointPool::append_state_hash(detail::StateHash& hash) const {
    hash.add(impl_->live_count);
    hash.add(static_cast<u64>(impl_->slots.size()));
    for (const auto& slot : impl_->slots) {
        hash.add(slot.generation); hash.add(slot.active);
        if (!slot.active)
            continue;
        const auto& r = slot.registration;
        hash.add(r.body_a); hash.add(r.body_b);
        for (const auto& edge : r.edges) {
            hash.add(edge.prev); hash.add(edge.next);
        }
        hash.add(r.type_index); hash.add(r.type); hash.add(r.collide_connected);
    }
    hash.sequence(impl_->free_list);
    const auto dense = [&](const auto& values) {
        hash.sequence(values.owner_slots);
        hash.add(static_cast<u64>(values.values.size()));
        for (const auto& joint : values.values) {
            hash.add(joint.local_anchor_a); hash.add(joint.local_anchor_b);
            if constexpr (requires { joint.local_axis_a; }) hash.add(joint.local_axis_a);
            if constexpr (requires { joint.reference_angle; }) hash.add(joint.reference_angle);
            if constexpr (requires { joint.lower_translation; }) {
                hash.add(joint.lower_translation); hash.add(joint.upper_translation);
                hash.add(joint.enable_limit);
            }
            if constexpr (requires { joint.length; }) hash.add(joint.length);
            if constexpr (requires { joint.rest_length; }) {
                hash.add(joint.rest_length); hash.add(joint.hertz);
                hash.add(joint.damping_ratio);
            }
        }
    };
    dense(impl_->pins); dense(impl_->sliders); dense(impl_->welds);
    dense(impl_->distances); dense(impl_->springs);
}

} // namespace sqp
