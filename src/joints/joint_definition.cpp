// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "sqphysics2d/joint.h"

#include "sqphysics2d/detail/data_storage.h"

namespace sqp {
namespace {

bool valid_common(const JointDef& def, const BodyPool& bodies) {
    if (!is_valid(def.body_b) || !bodies.is_valid(def.body_b) || !def.local_anchor_a.isfinite() ||
        !def.local_anchor_b.isfinite())
        return false;

    if (def.body_a == BodyHandle{})
        return def.local_anchor_a == Vec2{};

    return is_valid(def.body_a) && bodies.is_valid(def.body_a) && def.body_a != def.body_b;
}

} // namespace

bool validate_joint_def(const PinJointDef& def, const BodyPool& bodies) {
    return valid_common(def, bodies);
}

bool validate_joint_def(const SliderJointDef& def, const BodyPool& bodies) {
    return valid_common(def, bodies) && def.local_axis_a.isfinite() &&
           (def.local_axis_a.x != 0 || def.local_axis_a.y != 0) &&
           math::isfinite(def.lower_translation) && math::isfinite(def.upper_translation) &&
           def.lower_translation <= def.upper_translation;
}

bool validate_joint_def(const WeldJointDef& def, const BodyPool& bodies) {
    return valid_common(def, bodies);
}

bool validate_joint_def(const DistanceJointDef& def, const BodyPool& bodies) {
    return valid_common(def, bodies) && math::isfinite(def.length) && def.length >= 0;
}

bool validate_joint_def(const SpringJointDef& def, const BodyPool& bodies) {
    return valid_common(def, bodies) && math::isfinite(def.rest_length) && def.rest_length >= 0 &&
           math::isfinite(def.hertz) && def.hertz > 0 && math::isfinite(def.damping_ratio) &&
           def.damping_ratio >= 0;
}

} // namespace sqp
