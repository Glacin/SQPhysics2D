#include "scenario_runner.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <utility>

namespace sqp::benchmark {
namespace {

class Progress {
public:
    Progress() {
        show("preparing");
    }
    ~Progress() {
        std::fputc('\n', stderr);
    }

    void show(const std::string& message) {
        std::fprintf(stderr, "\r%s", message.c_str());
        if (previous_length_ > message.size())
            std::fprintf(stderr, "%*s", static_cast<int>(previous_length_ - message.size()), "");
        std::fflush(stderr);
        previous_length_ = message.size();
    }

    void show_tick(const char* phase, tick_t current, tick_t total) {
        show(std::string(phase) + " " + std::to_string(current) + "/" + std::to_string(total));
    }

private:
    std::size_t previous_length_ = 0;
};

#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
struct TickTotals {
    u64 intervals = 0;
    u64 active_islands = 0;
    u64 active_bodies = 0;
    u64 candidates = 0;
    u64 pairs = 0;
    u64 query_nanoseconds = 0;
    u64 update_nanoseconds = 0;
    u64 rebuild_nanoseconds = 0;
    u64 nodes_visited = 0;
    u64 proxy_reinsertions = 0;
    u64 query_calls = 0;

    void add(const WorldState::StageStats::IntervalEntry& entry) {
        ++intervals;
        active_islands += entry.active_islands;
        active_bodies += entry.active_bodies;
        candidates += entry.candidates;
        pairs += entry.pairs;
        query_nanoseconds += entry.query_nanoseconds;
        update_nanoseconds += entry.update_nanoseconds;
        rebuild_nanoseconds += entry.rebuild_nanoseconds;
        nodes_visited += entry.nodes_visited;
        proxy_reinsertions += entry.proxy_reinsertions;
        query_calls += entry.query_calls;
    }

    void add(const TickTotals& tick) {
        intervals += tick.intervals;
        active_islands += tick.active_islands;
        active_bodies += tick.active_bodies;
        candidates += tick.candidates;
        pairs += tick.pairs;
        query_nanoseconds += tick.query_nanoseconds;
        update_nanoseconds += tick.update_nanoseconds;
        rebuild_nanoseconds += tick.rebuild_nanoseconds;
        nodes_visited += tick.nodes_visited;
        proxy_reinsertions += tick.proxy_reinsertions;
        query_calls += tick.query_calls;
    }
};
struct QualityTotals {
    u64 height = 0;
    u64 reachable_nodes = 0;
    u64 leaf_count = 0;
    double area_ratio = 0;
    double sibling_overlap_ratio = 0;
};
#endif

} // namespace

double percentile_ms(std::vector<double> values, double percentile) {
    if (values.empty())
        return 0;
    std::sort(values.begin(), values.end());
    const auto rank = static_cast<std::size_t>(std::ceil(percentile * values.size()));
    return values[std::max<std::size_t>(1, rank) - 1];
}

Result run(Scenario& scenario, const Options& options) {
    if (options.measure_ticks == 0)
        throw std::invalid_argument("ticks must be positive");
    Progress progress;
    Context context;
    scenario.setup(context, options);
    if (options.proxy_extension_ratio)
        context.settings.broad_phase.proxy_extension_ratio = *options.proxy_extension_ratio;
    context.settings.allow_sleep = options.allow_sleep;
    if (options.toi_budget)
        context.settings.interval.max_truncations_per_tick = *options.toi_budget;
    context.settings.sanitize();
    context.state.settings = context.settings;
    std::vector<TickHash> tick_hashes;
    context.state.stats.record_body_freezes = options.freeze_stats;

    const auto step = [&](tick_t tick, const char* phase) {
        if (!context.world.step(config::TU_PER_TICK))
            throw std::runtime_error(std::string(phase) + " tick " + std::to_string(tick) +
                                     ": World::step failed");
        if (options.check_determinism) {
            const auto produced = context.world.last_step_tick_hashes();
            if (produced.size() != 1)
                throw std::runtime_error("expected one tick hash");
            tick_hashes.push_back(produced.front());
        }
    };
    const auto validate = [&](tick_t tick, const char* phase) {
        if (const auto error = scenario.validate(context); !error.empty())
            throw std::runtime_error(std::string(phase) + " tick " + std::to_string(tick) + ": " +
                                     error);
        if (const auto error = scenario.validate_sleep(context, tick); !error.empty())
            throw std::runtime_error(std::string(phase) + " tick " + std::to_string(tick) + ": " +
                                     error);
    };

    progress.show_tick("warmup", 0, options.warmup_ticks);
    for (tick_t tick = 1; tick <= options.warmup_ticks; ++tick) {
        step(tick, "warmup");
        validate(tick, "warmup");
        progress.show_tick("warmup", tick, options.warmup_ticks);
    }

    const tick_t measure_start = context.world.time_tu();
    auto body_freezes = std::move(context.state.stats.body_freezes);
    std::erase_if(body_freezes,
                  [measure_start](const auto& freeze) { return freeze.end <= measure_start; });
    context.state.stats = {};
    context.state.stats.record_toi_outcomes = options.toi_stats;
    context.state.stats.record_body_freezes = options.freeze_stats;
    context.state.stats.body_freezes = std::move(body_freezes);

    std::vector<double> times;
    times.reserve(static_cast<std::size_t>(options.measure_ticks));
    double total_seconds = 0;
    std::vector<real> half_layer_drift;
    std::vector<real> half_center_sag;
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
    TickTotals measured_totals;
    std::size_t previous_interval_count = 0;
    std::array<QualityTotals, 3> quality_totals{};
#endif
    progress.show_tick("ticks", 0, options.measure_ticks);
    for (tick_t tick = 1; tick <= options.measure_ticks; ++tick) {
        const auto begin = std::chrono::steady_clock::now();
        const bool stepped = context.world.step(config::TU_PER_TICK);
        const auto end = std::chrono::steady_clock::now();
        if (!stepped)
            throw std::runtime_error("measure tick " + std::to_string(tick) +
                                     ": World::step failed");
        if (options.check_determinism) {
            const auto produced = context.world.last_step_tick_hashes();
            if (produced.size() != 1)
                throw std::runtime_error("expected one tick hash");
            tick_hashes.push_back(produced.front());
        }
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
        TickTotals tick_totals;
        const auto& entries = context.state.stats.interval_entries;
        for (std::size_t i = previous_interval_count; i < entries.size(); ++i)
            tick_totals.add(entries[i]);
        previous_interval_count = entries.size();
        measured_totals.add(tick_totals);
        if (options.detail) {
            std::array<DynamicAabbTree::QualityStats, 3> quality;
            if (!context.state.broad_phase.tree.inspect(quality))
                throw std::runtime_error("broad phase tree integrity check failed");
            for (std::size_t type = 0; type < quality.size(); ++type) {
                quality_totals[type].height += quality[type].height;
                quality_totals[type].reachable_nodes += quality[type].reachable_nodes;
                quality_totals[type].leaf_count += quality[type].leaf_count;
                quality_totals[type].area_ratio += quality[type].area_ratio;
                quality_totals[type].sibling_overlap_ratio +=
                    quality[type].sibling_overlap_ratio;
            }
        }
#endif
        const double seconds = std::chrono::duration<double>(end - begin).count();
        total_seconds += seconds;
        times.push_back(seconds * 1000);
        validate(options.warmup_ticks + tick, "measure");
        if (tick == std::max<tick_t>(1, options.measure_ticks / 2)) {
            half_layer_drift = scenario.layer_horizontal_drift(context);
            half_center_sag = scenario.layer_center_sag(context);
        }
        progress.show_tick("ticks", tick, options.measure_ticks);
    }

    Result result;
    result.name = scenario.name();
    result.tick_hashes = std::move(tick_hashes);
    result.ticks_per_second =
        total_seconds > 0 ? static_cast<double>(options.measure_ticks) / total_seconds : 0;
    result.avg_ms_tick = total_seconds * 1000.0 / static_cast<double>(options.measure_ticks);
    result.simulation_percent =
        100.0 * result.ticks_per_second / static_cast<double>(config::TICKS_PER_SECOND);
    result.p50_tick_ms = percentile_ms(times, 0.50);
    result.p75_tick_ms = percentile_ms(times, 0.75);
    result.p95_tick_ms = percentile_ms(times, 0.95);
    result.max_tick_ms = *std::max_element(times.begin(), times.end());
    result.total_step_ms = total_seconds * 1000.0;
    const tick_t measure_end = context.world.time_tu();
    auto& freezes = context.state.stats.body_freezes;
    std::erase_if(freezes,
                  [measure_end](const auto& freeze) { return freeze.start >= measure_end; });
    for (auto& freeze : freezes) {
        freeze.start = std::max(freeze.start, measure_start);
        freeze.end = std::min(freeze.end, measure_end);
    }
    result.stage_stats = std::move(context.state.stats);
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
    const double tick_count = static_cast<double>(options.measure_ticks);
    result.tick_averages.intervals = measured_totals.intervals / tick_count;
    result.tick_averages.active_islands = measured_totals.active_islands / tick_count;
    result.tick_averages.active_bodies = measured_totals.active_bodies / tick_count;
    result.tick_averages.candidates = measured_totals.candidates / tick_count;
    result.tick_averages.pairs = measured_totals.pairs / tick_count;
    result.tick_averages.query_ms =
        static_cast<double>(measured_totals.query_nanoseconds) / 1'000'000.0 / tick_count;
    result.tick_averages.update_ms =
        static_cast<double>(measured_totals.update_nanoseconds) / 1'000'000.0 / tick_count;
    result.tick_averages.rebuild_ms =
        static_cast<double>(measured_totals.rebuild_nanoseconds) / 1'000'000.0 / tick_count;
    result.tick_averages.nodes_visited = measured_totals.nodes_visited / tick_count;
    result.tick_averages.proxy_reinsertions = measured_totals.proxy_reinsertions / tick_count;
    result.tick_averages.query_calls = measured_totals.query_calls / tick_count;
    result.tick_averages.nodes_per_query = measured_totals.query_calls
        ? static_cast<double>(measured_totals.nodes_visited) / measured_totals.query_calls : 0;
    if (options.detail)
        for (std::size_t type = 0; type < quality_totals.size(); ++type) {
            result.tree_quality[type].height =
                static_cast<u32>(quality_totals[type].height / options.measure_ticks);
            result.tree_quality[type].reachable_nodes =
                static_cast<u32>(quality_totals[type].reachable_nodes / options.measure_ticks);
            result.tree_quality[type].leaf_count =
                static_cast<u32>(quality_totals[type].leaf_count / options.measure_ticks);
            result.tree_quality[type].area_ratio = quality_totals[type].area_ratio / tick_count;
            result.tree_quality[type].sibling_overlap_ratio =
                quality_totals[type].sibling_overlap_ratio / tick_count;
        }
#endif
    result.half_layer_drift = std::move(half_layer_drift);
    result.final_layer_drift = scenario.layer_horizontal_drift(context);
    result.half_center_sag = std::move(half_center_sag);
    result.final_center_sag = scenario.layer_center_sag(context);
    return result;
}

} // namespace sqp::benchmark
