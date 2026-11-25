// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "world/world_stages.h"

#include "solver/contact_material.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace sqp {
namespace {

bool valid_material(const PhysicsMaterial& m) {
    return std::isfinite(m.contact_frequency_hz) && m.contact_frequency_hz >= real(0) &&
           std::isfinite(m.max_penetration_correction_speed) &&
           m.max_penetration_correction_speed >= real(0) &&
           std::isfinite(m.contact_damping_ratio) && m.contact_damping_ratio >= real(0) &&
           std::isfinite(m.base_restitution) && m.base_restitution >= real(0) &&
           m.base_restitution <= real(1) && std::isfinite(m.restitution_impact_threshold) &&
           m.restitution_impact_threshold > real(0) &&
           std::isfinite(m.restitution_impact_exponent) &&
           m.restitution_impact_exponent >= real(0) && std::isfinite(m.dynamic_friction) &&
           m.dynamic_friction >= real(0) && std::isfinite(m.static_friction) &&
           m.static_friction >= real(0) &&
           (m.restitution_flags & ~RESTITUTION_DISABLE_IMPACT_SCALING) == 0;
}

bool valid_manifold(const ContactManifold& manifold, bool allow_disabled) {
    if (manifold.count == 0)
        return allow_disabled;
    if (manifold.count < 1 || manifold.count > 2 || !std::isfinite(manifold.normal.x) ||
        !std::isfinite(manifold.normal.y))
        return false;
    const real length_squared = manifold.normal.length_squared();
    if (!std::isfinite(length_squared) || length_squared <= config::EPSILON * config::EPSILON)
        return false;
    for (i32 i = 0; i < manifold.count; ++i) {
        const auto& point = manifold.points[i];
        if (!std::isfinite(point.point.x) || !std::isfinite(point.point.y) ||
            !std::isfinite(point.depth))
            return false;
    }
    return true;
}

bool same_geometry(const ContactManifold& a, const ContactManifold& b) {
    if (a.count != b.count || a.normal.x != b.normal.x || a.normal.y != b.normal.y)
        return false;
    for (i32 i = 0; i < a.count; ++i)
        if (a.points[i].point.x != b.points[i].point.x ||
            a.points[i].point.y != b.points[i].point.y || a.points[i].depth != b.points[i].depth ||
            a.points[i].key != b.points[i].key)
            return false;
    return true;
}

bool valid_event_parameters(const PreSolveEvent& event) {
    return std::isfinite(event.combined_restitution) && event.combined_restitution >= real(0) &&
           event.combined_restitution <= real(1) &&
           std::isfinite(event.combined_dynamic_friction) &&
           event.combined_dynamic_friction >= real(0) &&
           std::isfinite(event.combined_static_friction) &&
           event.combined_static_friction >= real(0) &&
           event.combined_static_friction >= event.combined_dynamic_friction &&
           std::isfinite(event.combined_damping_ratio) && event.combined_damping_ratio >= real(0) &&
           std::isfinite(event.combined_contact_frequency) &&
           event.combined_contact_frequency >= real(0) &&
           std::isfinite(event.combined_max_penetration_correction_speed) &&
           event.combined_max_penetration_correction_speed >= real(0);
}

bool valid_layout(const WorldState::PairsData& data) {
    if (!data.pairs || !data.manifolds)
        return false;
    const std::size_t count = data.pairs->size();
    if (!data.toi_local_contacts.empty() && data.toi_local_contacts.size() != count)
        return false;
    if (count > std::numeric_limits<std::size_t>::max() / 2)
        return false;
    const auto& m = *data.manifolds;
    const std::size_t points = count * 2;
    return m.counts.size() == count && m.normals_x.size() == count && m.normals_y.size() == count &&
           m.keys.size() == points && m.points_x.size() == points && m.points_y.size() == points &&
           m.depths.size() == points;
}

bool load_body(const WorldState& state, BodyHandle handle, BodyState& body,
               BodyProperties& properties, BodyConfiguration& configuration) {
    return state.body_pool.is_valid(handle) && state.body_pool.get_state(handle, body) &&
           body.is_alive() && state.body_pool.get_properties(handle, properties) &&
           state.body_pool.get_configuration(handle, configuration);
}

real unlocked_inverse_mass(const BodyProperties& properties, const BodyConfiguration& config) {
    const auto flags = config.configs.override_flags;
    const bool locked_x =
        (flags & BodyConfigs::OverrideFlags::LockLinearX) != BodyConfigs::OverrideFlags(0);
    const bool locked_y =
        (flags & BodyConfigs::OverrideFlags::LockLinearY) != BodyConfigs::OverrideFlags(0);
    return locked_x && locked_y ? real(0) : properties.inv_mass;
}

real final_restitution(const WorldState& state, const ContactManifold& manifold,
                       const WorldState::ToiLocalContact* local,
                       const detail::CombinedContactMaterial& material, const BodyState& a,
                       const BodyState& b, const BodyProperties& properties_a,
                       const BodyProperties& properties_b, const BodyConfiguration& config_a,
                       const BodyConfiguration& config_b, const ShapeProperties& shape_a,
                       const ShapeProperties& shape_b) {

    const Vec2 normal = manifold.normal / std::sqrt(manifold.normal.length_squared());
    const Vec2 point = local && local->valid ? local->manifold.points[0].point
                                               : manifold.points[0].point;
    const Vec2 velocity_a = local && local->valid
                                ? a.linear_velocity + math::cross(
                                                          a.angular_velocity,
                                                          point - (a.position - local->origin))
                                : a.velocity_at_world_point(point);
    const Vec2 velocity_b = local && local->valid
                                ? b.linear_velocity + math::cross(
                                                          b.angular_velocity,
                                                          point - (b.position - local->origin))
                                : b.velocity_at_world_point(point);
    const real closing = -math::dot(velocity_b - velocity_a, normal);
    const real inverse_mass = unlocked_inverse_mass(properties_a, config_a) +
                              unlocked_inverse_mass(properties_b, config_b);
    return detail::contact_restitution(
        material, closing, inverse_mass,
        std::min(shape_a.restitution_impact_scale, shape_b.restitution_impact_scale),
        state.settings.contact.enable_dynamic_restitution);
}

} // namespace

bool run_pre_solve(WorldState& state) {
    auto& data = state.pairs;
    if (!valid_layout(data))
        return stage_error(state, "contact pair and manifold columns have invalid layout");

    const std::size_t count = data.pairs->size();
    if (data.original_manifolds_valid && data.original_manifolds.size() != count)
        return stage_error(state, "run_pre_solve: check failed (data.original_manifolds_valid && data.original_manifolds.size() != count)");
    auto& [manifolds, local_contacts, restitution, dynamic_friction, static_friction,
           damping, frequency, correction_speed] = data.pre_solve_scratch;
    manifolds = *data.manifolds;
    local_contacts = data.toi_local_contacts;
    for (auto* values : {&restitution, &dynamic_friction, &static_friction, &damping,
                         &frequency, &correction_speed})
        values->resize(count);

    for (std::size_t i = 0; i < count; ++i) {
        const Pair pair = (*data.pairs)[i];
        if (pair.body_a == pair.body_b || pair.shape_a == pair.shape_b)
            return stage_error(state, "run_pre_solve: check failed (pair.body_a == pair.body_b || pair.shape_a == pair.shape_b)");

        BodyState body_a, body_b;
        BodyProperties properties_a, properties_b;
        BodyConfiguration config_a, config_b;
        BodyHandle owner_a, owner_b;
        ShapeProperties shape_a, shape_b;
        PhysicsMaterial material_a, material_b;

        // Careful Validation
        if (!load_body(state, pair.body_a, body_a, properties_a, config_a) ||
            !load_body(state, pair.body_b, body_b, properties_b, config_b) ||
            !state.shape_pool.is_valid(pair.shape_a) || !state.shape_pool.is_valid(pair.shape_b) ||
            !state.body_pool.get_shape_owner(pair.shape_a, owner_a) ||
            !state.body_pool.get_shape_owner(pair.shape_b, owner_b) || owner_a != pair.body_a ||
            owner_b != pair.body_b || !state.shape_pool.get_properties(pair.shape_a, shape_a) ||
            !state.shape_pool.get_properties(pair.shape_b, shape_b) ||
            !state.shape_pool.get_material(pair.shape_a, material_a) ||
            !state.shape_pool.get_material(pair.shape_b, material_b) ||
            !valid_material(material_a) || !valid_material(material_b) ||
            !std::isfinite(shape_a.restitution_impact_scale) ||
            !std::isfinite(shape_b.restitution_impact_scale) ||
            shape_a.restitution_impact_scale <= real(0) ||
            shape_b.restitution_impact_scale <= real(0))
            return stage_error(state, "run_pre_solve: contact pair has invalid body, shape, or shape owner");

        ContactManifold manifold = manifolds.at(i);
        if (!valid_manifold(manifold, true))
            return stage_error(state, "run_pre_solve: check failed (!valid_manifold(manifold, true))");
        const auto combined = detail::combine_contact_materials(material_a, material_b);
        PreSolveEvent result;
        result.contact_pair = pair;
        result.manifold = manifold;
        result.combined_restitution =
            manifold.count == 0
                ? real(0)
                : final_restitution(state, manifold,
                                    local_contacts.empty() ? nullptr : &local_contacts[i],
                                    combined, body_a, body_b, properties_a,
                                    properties_b, config_a, config_b, shape_a, shape_b);
        result.combined_dynamic_friction = combined.dynamic_friction;
        result.combined_static_friction = combined.static_friction;
        result.combined_damping_ratio = combined.contact_damping_ratio;
        result.combined_contact_frequency = combined.contact_frequency_hz;
        result.combined_max_penetration_correction_speed =
            combined.max_penetration_correction_speed;

        const auto flag = BodyConfigs::EventFlags::PreSolve;
        bool callback_used = false;
        if (state.callbacks.pre_solve &&
            (((config_a.configs.event_flags & flag) != BodyConfigs::EventFlags(0)) ||
             ((config_b.configs.event_flags & flag) != BodyConfigs::EventFlags(0)))) {
            PhysicsEvent event{PhysicsEventTypes::PreSolve, result};
            state.callbacks.pre_solve(event, state.callbacks.pre_data);
            const auto* updated = std::get_if<PreSolveEvent>(&event.event);
            if (event.type != PhysicsEventTypes::PreSolve || !updated ||
                updated->contact_pair != pair)
                return stage_error(state, "run_pre_solve: check failed (event.type != PhysicsEventTypes::PreSolve || !updated || updated->contact_pair != pair)");
            result = *updated;
            callback_used = true;
        }
        bool validate_result = callback_used;
#ifdef SQP_DEBUG_VALIDATION
        validate_result = true;
#endif
        if (validate_result &&
            (!valid_manifold(result.manifold, true) || !valid_event_parameters(result)))
            return stage_error(state, "run_pre_solve: check failed (validate_result && (!valid_manifold(result.manifold, true) || !valid_event_parameters(result)))");

        if (!local_contacts.empty() && local_contacts[i].valid &&
            !same_geometry(manifold, result.manifold)) {
            local_contacts[i].manifold = result.manifold;
            for (i32 slot = 0; slot < result.manifold.count; ++slot)
                local_contacts[i].manifold.points[slot].point -= local_contacts[i].origin;
        }

        const ContactManifold& output = result.manifold;
        manifolds.counts[i] = static_cast<u8>(output.count);
        manifolds.normals_x[i] = output.normal.x;
        manifolds.normals_y[i] = output.normal.y;
        for (std::size_t slot = 0; slot < 2; ++slot) {
            const std::size_t index = i * 2 + slot;
            manifolds.keys[index] = output.points[slot].key;
            manifolds.points_x[index] = output.points[slot].point.x;
            manifolds.points_y[index] = output.points[slot].point.y;
            manifolds.depths[index] = output.points[slot].depth;
        }
        restitution[i] = result.combined_restitution;
        dynamic_friction[i] = result.combined_dynamic_friction;
        static_friction[i] = result.combined_static_friction;
        damping[i] = result.combined_damping_ratio;
        frequency[i] = result.combined_contact_frequency;
        correction_speed[i] = result.combined_max_penetration_correction_speed;
    }

    if (!data.original_manifolds_valid) {
        data.original_manifolds = *data.manifolds;
        data.original_manifolds_valid = true;
    }
    std::swap(*data.manifolds, manifolds);
    data.toi_local_contacts.swap(local_contacts);
    data.combined_restitutions.swap(restitution);
    data.combined_dynamic_frictions.swap(dynamic_friction);
    data.combined_static_frictions.swap(static_friction);
    data.combined_damping_ratios.swap(damping);
    data.combined_contact_frequencies.swap(frequency);
    data.combined_max_penetration_correction_speeds.swap(correction_speed);
    return true;
}

} // namespace sqp
