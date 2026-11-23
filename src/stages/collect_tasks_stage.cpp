// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "world/world_stages.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <variant>

namespace sqp {
namespace {

template <class F> bool each_body(const Island& island, F&& fn) {
    if (const auto* small = std::get_if<Island::SmallIsland>(&island.bodies)) {
        for (BodyHandle body : small->data) {
            if (!is_valid(body))
                break;
            if (!fn(body))
                return false;
        }
    } else {
        for (BodyHandle body : std::get<std::vector<BodyHandle>>(island.bodies))
            if (!fn(body))
                return false;
    }
    return true;
}

bool uses_ccd_for_tick(const WorldState& state, BodyHandle body,
                       const BodyConfiguration& configuration, real speed, real turn,
                       std::vector<ShapeHandle>& shapes) {
    const auto mode = configuration.configs.ccd_override;
    if (mode == BodyConfigs::CcdOverride::Disabled ||
        (!state.settings.ccd.enable_circle && !state.settings.ccd.enable_polygon))
        return false;
    const real tick_seconds = math::tu_to_seconds<real>(config::TU_PER_TICK);
    const real threshold = state.settings.ccd.linear_start_threshold_m;
    const bool linear_cast = speed * tick_seconds > threshold;
    if (mode == BodyConfigs::CcdOverride::Auto && !linear_cast && turn == real(0))
        return false;
    if (!state.body_pool.get_attached_shapes(body, shapes))
        return false;
    for (ShapeHandle shape : shapes) {
        if ((shape.kind == ShapeKind::Circle && !state.settings.ccd.enable_circle) ||
            (shape.kind == ShapeKind::Polygon && !state.settings.ccd.enable_polygon))
            continue;
        if (mode == BodyConfigs::CcdOverride::Shapecast || linear_cast ||
            turn * tick_seconds * state.shape_pool.sweep_radius(shape) > threshold)
            return true;
    }
    return false;
}

bool period_and_sleep(WorldState& state, IslandId id, tick_t elapsed, tick_t& period,
                      bool& sleeping, bool& uses_ccd) {
    const Island* island = state.island_pool.get(id);
    if (!island)
        return stage_error(state, "period_and_sleep: check failed (!island)");
    const auto& limits = state.settings.interval;
    period = config::TU_PER_TICK;
    bool quiet = true;
    bool no_sleep = false;
    uses_ccd = false;
    std::vector<ShapeHandle> shapes;
    const bool valid = each_body(*island, [&](BodyHandle body) {
        BodyState value;
        BodyConfiguration config;
        if (!state.body_pool.get_state(body, value) ||
            !state.body_pool.get_configuration(body, config) ||
            config.configs.register_type == BodyConfigs::RegisterType::ImmutableStatic ||
            !value.linear_velocity.isfinite() || !std::isfinite(value.angular_velocity))
            return stage_error(state, "period_and_sleep: check failed (!state.body_pool.get_state(body, value) || !state.body_pool.get_configuration(body, config) || config.configs.register_type == BodyConfigs::RegisterType::ImmutableStatic || !value.linear_velocity.isfinite() || !std::isfinite(value.angular_velocity))");
#ifdef SQP_DEBUG_VALIDATION
        if (state.body_pool.island_of(body) != id)
            return stage_error(state, "period_and_sleep: check failed (state.body_pool.island_of(body) != id)");
#endif
        const real speed = value.linear_velocity.length();
        const real turn = math::abs(value.angular_velocity);
        const auto candidate = [&](real cap, real velocity) {
            if (!(velocity > real(0)) || !(cap > real(0)))
                return config::TU_PER_TICK;
            const long double raw = static_cast<long double>(cap) *
                                    static_cast<long double>(config::TU_PER_SECOND) /
                                    static_cast<long double>(velocity);
            if (raw >= static_cast<long double>(config::TU_PER_TICK))
                return config::TU_PER_TICK;
            return std::max(limits.min_period_tu,
                            std::bit_floor(std::max(tick_t(1), static_cast<tick_t>(raw))));
        };
        const bool body_uses_ccd = uses_ccd_for_tick(state, body, config, speed, turn, shapes);
        uses_ccd |= body_uses_ccd;
        if (!body_uses_ccd) {
            period = std::min(period, candidate(limits.max_linear_displacement_m, speed));
            period = std::min(period, candidate(limits.max_angular_displacement_rad, turn));
        }
        if (config.configs.register_type == BodyConfigs::RegisterType::Kinematic) {
            quiet &= speed == real(0) && turn == real(0);
        } else {
            quiet &= speed <= state.settings.sleep_linear_velocity_threshold_m_s &&
                     turn <= state.settings.sleep_angular_velocity_threshold_rad_s;
        }
        no_sleep |= math::has_any_flags(config.configs.override_flags,
                                         BodyConfigs::OverrideFlags::NoSleep);
        return true;
    });
    if (!valid)
        return stage_error(state, "period_and_sleep: check failed (!valid)");
    sleeping = false;
    if (state.island_pool.frozen_until(id) > state.schedule.time)
        return state.island_pool.set_sleep_time(id, 0) ||
               stage_error(state, "period_and_sleep: failed to reset frozen island sleep time");
    if (!state.settings.allow_sleep || !quiet || no_sleep)
        return state.island_pool.set_sleep_time(id, 0) ||
               stage_error(state, "period_and_sleep: failed to reset island sleep time");
    const tick_t previous = state.island_pool.sleep_time(id);
    const tick_t total = elapsed > std::numeric_limits<tick_t>::max() - previous
                             ? std::numeric_limits<tick_t>::max()
                             : previous + elapsed;
    if (!state.island_pool.set_sleep_time(id, total))
        return stage_error(state, "period_and_sleep: check failed (!state.island_pool.set_sleep_time(id, total))");
    sleeping = total >= state.settings.sleep_time_threshold_tu;
    return true;
}

} // namespace

bool unschedule_island(WorldState& state, IslandId id) {
    u8 bucket;
    std::size_t position;
    if (!state.island_pool.get_schedule_location(id, bucket, position))
        return stage_error(state, "unschedule_island: check failed (!state.island_pool.get_schedule_location(id, bucket, position))");
    if (bucket == 255)
        return true;
    if (bucket >= state.schedule.period_islands.size())
        return stage_error(state, "unschedule_island: check failed (bucket >= state.schedule.period_islands.size())");
    auto& entries = state.schedule.period_islands[bucket];
    if (position >= entries.size() || entries[position] != id)
        return stage_error(state, "unschedule_island: check failed (position >= entries.size() || entries[position] != id)");
    if (position + 1 != entries.size()) {
        const IslandId moved = entries.back();
        entries[position] = moved;
        if (!state.island_pool.set_schedule_location(moved, bucket, position))
            return stage_error(state, "unschedule_island: check failed (!state.island_pool.set_schedule_location(moved, bucket, position))");
        state.schedule.entries[moved.index].position = position;
    }
    entries.pop_back();
    return state.island_pool.set_schedule_location(id, 255, 0) ||
           stage_error(state, "unschedule_island: failed to clear schedule location");
}

bool park_island(WorldState& state, IslandId id) {
    const Island* island = state.island_pool.get(id);
    if (!island || !unschedule_island(state, id))
        return stage_error(state, "park_island: check failed (!island || !unschedule_island(state, id))");
    if (!state.island_pool.set_sleep_time(id, state.settings.sleep_time_threshold_tu))
        return stage_error(state, "park_island: check failed (!state.island_pool.set_sleep_time(id, state.settings.sleep_time_threshold_tu))");
    if (!each_body(*island, [&](BodyHandle body) {
            BodyState value;
            if (!state.body_pool.get_state(body, value))
                return stage_error(state, "park_island: check failed (!state.body_pool.get_state(body, value))");
            value.linear_velocity = {};
            value.angular_velocity = 0;
            if (!state.body_pool.set_state(body, value))
                return stage_error(state, "park_island: check failed (!state.body_pool.set_state(body, value))");
            const auto dense = state.body_pool.dense_index(body);
            if (!dense)
                return stage_error(state, "park_island: check failed (!dense)");
            ++state.body_pool.motion_version_column(body.kind())[*dense];
            return true;
        }))
        return stage_error(state, "park_island: check failed (!dense)");
    auto& bucket = state.schedule.period_islands[0];
    const std::size_t position = bucket.size();
    bucket.push_back(id);
    if (id.index >= state.schedule.entries.size())
        state.schedule.entries.resize(std::size_t(id.index) + 1);
    state.schedule.entries[id.index] = {id, 0, 0, position};
    return state.island_pool.set_schedule_location(id, 0, position) ||
           stage_error(state, "park_island: failed to set sleeping schedule location");
}

bool reschedule_island(WorldState& state, IslandId id, tick_t elapsed) {
    const Island* island = state.island_pool.get(id);
    if (!island)
        return stage_error(state, "reschedule_island: check failed (!island)");
    tick_t period;
    bool sleeping;
    bool uses_ccd;
    if (!period_and_sleep(state, id, elapsed, period, sleeping, uses_ccd) || period == 0)
        return stage_error(state, "reschedule_island: check failed (!period_and_sleep(state, id, elapsed, period, sleeping, uses_ccd) || period == 0)");
    if (sleeping)
        return park_island(state, id);
    const unsigned index = std::countr_zero(period) + 1;
    if (index >= state.schedule.period_islands.size() || !unschedule_island(state, id))
        return stage_error(state, "reschedule_island: check failed (index >= state.schedule.period_islands.size() || !unschedule_island(state, id))");
    auto& bucket = state.schedule.period_islands[index];
    const std::size_t position = bucket.size();
    bucket.push_back(id);
    if (id.index >= state.schedule.entries.size())
        state.schedule.entries.resize(std::size_t(id.index) + 1);
    state.schedule.entries[id.index] = {id, period, 0, position, uses_ccd};
    return state.island_pool.set_schedule_location(id, static_cast<u8>(index), position) ||
           stage_error(state, "reschedule_island: failed to store schedule location");
}

bool wake_island(WorldState& state, IslandId id) {
    u8 bucket;
    std::size_t position;
    if (!state.island_pool.get_schedule_location(id, bucket, position))
        return stage_error(state, "wake_island: check failed (!state.island_pool.get_schedule_location(id, bucket, position))");
    if (bucket == 0) {
        const Island* island = state.island_pool.get(id);
        if (!island || !state.island_pool.set_time(id, std::max(state.schedule.time,
                                                                island->time)))
            return stage_error(state, "wake_island: check failed (!island || !state.island_pool.set_time(id, std::max(state.schedule.time, island->time)))");
    }
    if (!state.island_pool.set_sleep_time(id, 0))
        return stage_error(state, "wake_island: failed to reset sleep time");
    return reschedule_island(state, id);
}

bool run_collect_tasks(WorldState& state) {
    if (state.interval.microstep)
        return true;
    if (!run_rebuild_islands(state, false))
        return stage_error(state, "run_collect_tasks: check failed (!run_rebuild_islands(state, false))");
    auto& schedule = state.schedule;
    for (const auto& processed : state.solve.processed_islands) {
        if (processed.end < processed.start)
            return stage_error(state, "processed island interval ends before it starts");
        if (!state.island_pool.get(processed.id))
            continue;
        if (!reschedule_island(state, processed.id, processed.end - processed.start))
            return stage_error(state, "run_collect_tasks: check failed (!reschedule_island(state, processed.id, processed.end - processed.start))");
    }
    state.solve.processed_islands.clear();
    tick_t next = std::numeric_limits<tick_t>::max();
    for (std::size_t index = 1; index < schedule.period_islands.size(); ++index) {
        if (schedule.period_islands[index].empty())
            continue;
        const tick_t period = tick_t(1) << (index - 1);
        const tick_t remainder = schedule.time % period;
        const tick_t base = schedule.time;
        const tick_t offset = remainder == 0 ? period : period - remainder;
        if (offset > std::numeric_limits<tick_t>::max() - base)
            return stage_error(state, "run_collect_tasks: check failed (offset > std::numeric_limits<tick_t>::max() - base)");
        next = std::min(next, base + offset);
    }
    // A frozen island is due at its resume boundary even if its usual period is shorter.
    for (std::size_t index = 1; index < schedule.period_islands.size(); ++index)
        for (IslandId id : schedule.period_islands[index])
            if (const tick_t until = state.island_pool.frozen_until(id);
                until > schedule.time)
                next = std::min(next, until);
    state.interval.islands->clear();
    if (next == std::numeric_limits<tick_t>::max()) {
        schedule.ordinary_target = next;
        state.interval.tick = false;
        schedule.prepared = true;
        return true;
    }
    for (std::size_t index = 1; index < schedule.period_islands.size(); ++index) {
        const tick_t period = tick_t(1) << (index - 1);
        for (IslandId id : schedule.period_islands[index])
            if (const Island* island = state.island_pool.get(id);
                island && island->time <= next &&
                (next % period == 0 || schedule.entries[id.index].uses_ccd))
                state.interval.islands->push_back(id);
    }
    schedule.ordinary_target = next;
    state.interval.target_end = next;
    state.interval.tick = next % config::TU_PER_TICK == 0;
    schedule.prepared = false;
    return true;
}

} // namespace sqp
