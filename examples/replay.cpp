// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include <sqphysics2d/build_config.h>
#include <sqphysics2d/world.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

constexpr std::size_t tick_count = 8;
constexpr auto tick = sqp::config::TU_PER_TICK;
using Replay = std::array<sqp::TickHash, tick_count>;

void step(sqp::World& world, sqp::tick_t duration) {
    if (!world.step(duration))
        throw std::runtime_error("Simulation step failed: " + std::string(world.last_error_message()));
}

sqp::CommandResult take_applied(sqp::World& world, sqp::CommandId id, std::string_view action) {
    sqp::CommandResult result;
    if (!world.take_command_result(id, result) || result.status != sqp::CommandStatus::Applied)
        throw std::runtime_error(std::string(action) + " command was not applied");
    return result;
}

Replay run(sqp::World& world, sqp::BodyHandle body) {
    const auto impulse = world.apply_impulse(world.time_tu() + 2 * tick, body, {1, 0});
    step(world, tick_count * tick);
    take_applied(world, impulse, "Impulse");

    const auto hashes = world.last_step_tick_hashes();
    if (hashes.size() != tick_count)
        throw std::runtime_error("Expected one hash per tick");

    // The span belongs to the world and is invalidated on restore.
    Replay replay;
    std::copy(hashes.begin(), hashes.end(), replay.begin());
    return replay;
}

} // namespace

int main() {
    try {
        sqp::World world;
        const auto create_body = world.create_rigid(0);
        step(world, 0);
        const auto body = take_applied(world, create_body, "Body creation").body;

        sqp::CircleShapeDef circle;
        circle.radius = 0.5;
        const auto create_shape = world.create_circle_shape(0, circle, body);
        step(world, 0);
        take_applied(world, create_shape, "Shape creation");

        sqp::WorldSnapshot snapshot;
        if (!world.capture_snapshot(snapshot))
            throw std::runtime_error("Snapshot capture failed");
        const auto original = run(world, body);
        if (!world.restore_snapshot(snapshot))
            throw std::runtime_error("Snapshot restore failed");
        const auto replay = run(world, body);

        for (std::size_t i = 0; i < tick_count; ++i) {
            if (original[i].tick != replay[i].tick || original[i].hash != replay[i].hash) {
                std::cerr << "Replay diverged at tick " << original[i].tick << ": expected "
                          << std::hex << original[i].hash << ", got " << replay[i].hash
                          << std::dec << " at tick " << replay[i].tick << '\n';
                return 1;
            }
        }

        std::cout << "Replay matched " << tick_count << " ticks.\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
