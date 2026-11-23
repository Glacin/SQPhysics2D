#include "scenario_runner.h"
#include "scenarios.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>

#ifdef _WIN32
#include <io.h>
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <unistd.h>
#endif

namespace {

bool parse_number(const char* text, sqp::u64& value) {
    if (!text || *text == '-' || *text == '\0')
        return false;
    errno = 0;
    char* end = nullptr;
    const auto parsed = std::strtoull(text, &end, 0);
    if (errno == ERANGE || end == text || *end != '\0')
        return false;
    value = static_cast<sqp::u64>(parsed);
    return true;
}

bool parse_ratio(const char* text, sqp::real& value) {
    if (!text || *text == '\0')
        return false;
    errno = 0;
    char* end = nullptr;
    const double parsed = std::strtod(text, &end);
    if (errno == ERANGE || end == text || *end != '\0' || !std::isfinite(parsed) ||
        parsed < 0 || parsed > 1)
        return false;
    value = static_cast<sqp::real>(parsed);
    return true;
}

void usage() {
    std::fprintf(stderr, "usage: sqphysics2d_benchmark SCENARIO [--bodies N] [--seed N] "
                         "[--warmup N] [--ticks N] [--sleep on|off] "
                         "[--initial-sleep on|off] [--awake-bodies N] "
                         "[--proxy-extension-ratio 0..1] [--detail] [--check-determinism]\n"
                         "scenarios: box_free_motion, box_floor_sleep, aabb_box_pyramid, "
                         "rigid_box_pyramid, circle_box_bounce, convex_polygon_bounce, "
                         "mixed_shape_bounce\n"
                         "--bodies defaults to 50 for pyramid bottom rows and 1000 otherwise.\n");
}

#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
bool color_output() {
#ifdef _WIN32
    if (!_isatty(_fileno(stdout)))
        return false;
    const HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
    DWORD mode = 0;
    return handle != INVALID_HANDLE_VALUE && GetConsoleMode(handle, &mode) &&
           SetConsoleMode(handle, mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
#else
    return isatty(fileno(stdout)) != 0;
#endif
}
#endif

void print_detail(const sqp::benchmark::Result& result) {
    constexpr std::array<const char*, sqp::WorldState::StageStats::stage_count> names = {
        "process_command", "collect_tasks", "broad_phase", "narrow_phase", "process_toi",
        "update_pairs", "pre_solve", "solve", "post_solve", "advance_interval"};
    std::array<std::size_t, names.size()> order{};
    for (std::size_t i = 0; i < order.size(); ++i)
        order[i] = i;
    const auto& stages = result.stage_stats.stages;
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return stages[a].nanoseconds > stages[b].nanoseconds;
    });
    const bool color = color_output();
#endif
    std::printf("\nintervals=%llu\n", static_cast<unsigned long long>(result.stage_stats.intervals));
    std::printf("%-20s %16s %12s %12s\n", "stage", "total_ms", "step_percent", "calls");
    for (std::size_t rank = 0; rank < order.size(); ++rank) {
        const std::size_t i = order[rank];
        const auto calls = static_cast<unsigned long long>(stages[i].calls);
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
        const double ms = static_cast<double>(stages[i].nanoseconds) / 1'000'000.0;
        const double percent = result.total_step_ms > 0 ? 100.0 * ms / result.total_step_ms : 0;
        const char* ink = color && stages[i].nanoseconds > 0
                              ? (rank == 0 ? "\x1b[31m" : rank < 3 ? "\x1b[33m" : "")
                              : "";
        std::printf("%s%-20s %16.6f %11.2f%% %12llu%s\n", ink, names[i], ms, percent,
                    calls, *ink ? "\x1b[0m" : "");
#else
        std::printf("%-20s %16s %12s %12llu\n", names[i], "n/a", "n/a", calls);
#endif
    }
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
    const auto& average = result.tick_averages;
    std::printf("\n%-16s %16s %16s %16s %16s %16s\n", "tick_average", "intervals",
                "active_islands", "active_bodies", "candidates", "pairs");
    std::printf("%-16s %16.2f %16.2f %16.2f %16.2f %16.2f\n", "per_tick",
                average.intervals, average.active_islands, average.active_bodies,
                average.candidates, average.pairs);
    std::printf("\n%-16s %16s %16s %16s %20s\n", "tick_average", "query_ms",
                "update_ms", "nodes_visited", "proxy_reinsertions");
    std::printf("%-16s %16.6f %16.6f %16.2f %20.2f\n", "per_tick", average.query_ms,
                average.update_ms, average.nodes_visited, average.proxy_reinsertions);
    std::printf("%-16s %16s %16s\n", "query_load", "queries/tick", "nodes/query");
    std::printf("%-16s %16.2f %16.2f\n", "per_tick", average.query_calls,
                average.nodes_per_query);
    std::printf("%-16s %16.6f\n", "rebuild_ms/tick", average.rebuild_ms);
    std::printf("\n%-12s %8s %12s %12s %14s %16s\n", "tree", "height",
                "nodes", "proxies", "area_ratio", "sibling_overlap");
    constexpr std::array<const char*, 3> types = {"static", "kinematic", "dynamic"};
    for (std::size_t type = 0; type < types.size(); ++type) {
        const auto& q = result.tree_quality[type];
        std::printf("%-12s %8u %12u %12u %14.4f %16.4f\n", types[type], q.height,
                    q.reachable_nodes, q.leaf_count, q.area_ratio,
                    q.sibling_overlap_ratio);
    }
#endif
}

} // namespace

int main(int argc, char** argv) {
    using namespace sqp;
    using namespace sqp::benchmark;
    if (argc < 2) {
        usage();
        return 1;
    }
    auto scenario = make_scenario(argv[1]);
    if (!scenario) {
        std::fprintf(stderr, "unknown scenario: %s\n", argv[1]);
        usage();
        return 1;
    }
    Options options;
    const std::string name = scenario->name();
    const bool pyramid = name == "aabb_box_pyramid" || name == "rigid_box_pyramid";
    if (pyramid)
        options.body_count = 50;
    options.allow_sleep = scenario->default_allow_sleep();

    for (int i = 2; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "--detail") {
            options.detail = true;
            continue;
        }
        if (argument == "--check-determinism") {
            options.check_determinism = true;
            continue;
        }
        if (argument != "--bodies" && argument != "--seed" && argument != "--warmup" &&
            argument != "--ticks" && argument != "--sleep" && argument != "--initial-sleep" &&
            argument != "--awake-bodies" && argument != "--proxy-extension-ratio") {
            std::fprintf(stderr, "unknown option: %s\n", argument.c_str());
            return 1;
        }
        if (i + 1 >= argc) {
            std::fprintf(stderr, "missing value for %s\n", argument.c_str());
            return 1;
        }
        const char* value_text = argv[++i];
        if (argument == "--proxy-extension-ratio") {
            real ratio = 0;
            if (!parse_ratio(value_text, ratio)) {
                std::fprintf(stderr, "proxy-extension-ratio must be between 0 and 1\n");
                return 1;
            }
            options.proxy_extension_ratio = ratio;
            continue;
        }
        if (argument == "--sleep" || argument == "--initial-sleep") {
            const std::string value = value_text;
            if (value != "on" && value != "off") {
                std::fprintf(stderr, "%s must be on or off\n", argument.c_str());
                return 1;
            }
            if (argument == "--sleep")
                options.allow_sleep = value == "on";
            else
                options.initial_sleep = value == "on";
            continue;
        }
        u64 value = 0;
        if (!parse_number(value_text, value)) {
            std::fprintf(stderr, "invalid numeric value for %s\n", argument.c_str());
            return 1;
        }
        if (argument == "--bodies" || argument == "--seed" || argument == "--awake-bodies") {
            if (value > std::numeric_limits<u32>::max()) {
                std::fprintf(stderr, "%s exceeds u32 range\n", argument.c_str());
                return 1;
            }
            if (argument == "--bodies")
                options.body_count = static_cast<u32>(value);
            else if (argument == "--seed")
                options.seed = static_cast<u32>(value);
            else
                options.awake_bodies = static_cast<u32>(value);
        } else if (argument == "--warmup") {
            options.warmup_ticks = value;
        } else if (argument == "--ticks") {
            options.measure_ticks = value;
        }
    }
    if (options.body_count == 0) {
        std::fprintf(stderr, "bodies must be positive\n");
        return 1;
    }
    if (name != "box_floor_sleep" && (options.initial_sleep || options.awake_bodies != 0)) {
        std::fprintf(stderr, "initial-sleep and awake-bodies require box_floor_sleep\n");
        return 1;
    }

    try {
        const Result result = run(*scenario, options);
        if (options.check_determinism) {
            auto replay = make_scenario(name);
            const Result replay_result = run(*replay, options);
            if (const auto mismatch = first_hash_mismatch(result.tick_hashes,
                                                           replay_result.tick_hashes)) {
                const std::size_t i = *mismatch;
                if (i >= result.tick_hashes.size() || i >= replay_result.tick_hashes.size())
                    throw std::runtime_error("determinism: tick count differs");
                const auto a = result.tick_hashes[i];
                const auto b = replay_result.tick_hashes[i];
                char message[180];
                std::snprintf(message, sizeof(message),
                              "determinism mismatch at tick %llu: %016llx vs %016llx",
                              static_cast<unsigned long long>(a.tick),
                              static_cast<unsigned long long>(a.hash),
                              static_cast<unsigned long long>(b.hash));
                throw std::runtime_error(message);
            }
            std::printf("determinism=pass ticks=%zu\n", result.tick_hashes.size());
        }
        std::printf("scenario=%s\nsimulation=%.2f%%\nticks_per_second=%.2f\n"
                    "avg_ms_tick=%.6f\np50_tick_ms=%.6f\np75_tick_ms=%.6f\np95_tick_ms=%.6f\n"
                    "max_tick_ms=%.6f\n",
                    result.name.c_str(), result.simulation_percent, result.ticks_per_second,
                    result.avg_ms_tick, result.p50_tick_ms, result.p75_tick_ms, result.p95_tick_ms,
                    result.max_tick_ms);
        if (options.detail) {
            for (std::size_t row = 0; row < result.final_layer_drift.size(); ++row)
                std::printf("layer_%zu_horizontal_drift_half=%.6f final=%.6f\n", row,
                            static_cast<double>(result.half_layer_drift[row]),
                            static_cast<double>(result.final_layer_drift[row]));
            for (std::size_t row = 0; row < result.final_center_sag.size(); ++row)
                std::printf("layer_%zu_center_sag_half=%.6f final=%.6f\n", row,
                            static_cast<double>(result.half_center_sag[row]),
                            static_cast<double>(result.final_center_sag[row]));
            print_detail(result);
        }
        return 0;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 2;
    }
}
