#pragma once

#include <algorithm>
#include <string>
#include <optional>
#include <span>
#include <vector>

#include "sqphysics2d/detail/core_world.h"

namespace sqp::benchmark {

struct Options {
    u32 body_count = 1000;
    u32 seed = 0x05162020u;
    tick_t warmup_ticks = 120;
    tick_t measure_ticks = 300;
    bool allow_sleep = true;
    bool initial_sleep = false;
    u32 awake_bodies = 0;
    bool detail = false;
    bool check_determinism = false;
    bool freeze_stats = false;
    bool toi_stats = false;
    std::optional<real> proxy_extension_ratio;
    std::optional<u32> toi_budget;
};

struct Context {
    World world;
    WorldState& state = world.state();
    PhysicsSettings settings;
};

struct PositionStats {
    bool available = false;
    real position_delta = 0;
    real position_limit = 0;
};

struct SleepStats {
    u32 awake = 0;
    u32 sleeping = 0;
};

struct Result {
    std::string name;
    std::vector<TickHash> tick_hashes;
    double simulation_percent = 0;
    double ticks_per_second = 0;
    double avg_ms_tick = 0;
    double p50_tick_ms = 0;
    double p75_tick_ms = 0;
    double p95_tick_ms = 0;
    double max_tick_ms = 0;
    double total_step_ms = 0;
    WorldState::StageStats stage_stats;
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
    struct TickAverages {
        double intervals = 0;
        double active_islands = 0;
        double active_bodies = 0;
        double candidates = 0;
        double pairs = 0;
        double query_ms = 0;
        double update_ms = 0;
        double rebuild_ms = 0;
        double nodes_visited = 0;
        double proxy_reinsertions = 0;
        double query_calls = 0;
        double nodes_per_query = 0;
    } tick_averages;
    std::array<DynamicAabbTree::QualityStats, 3> tree_quality{};
#endif
    std::vector<real> half_layer_drift, final_layer_drift;
    std::vector<real> half_center_sag, final_center_sag;
};

class Scenario {
public:
    virtual ~Scenario() = default;
    virtual const char* name() const = 0;
    virtual void setup(Context&, const Options&) = 0;
    virtual std::string validate(const Context&) const = 0;
    virtual std::vector<real> layer_horizontal_drift(const Context&) const {
        return {};
    }
    virtual std::vector<real> layer_center_sag(const Context&) const {
        return {};
    }
    virtual std::string validate_sleep(const Context&, tick_t) const {
        return {};
    }
    virtual bool default_allow_sleep() const {
        return true;
    }
};

inline std::optional<std::size_t> first_hash_mismatch(std::span<const TickHash> a,
                                                       std::span<const TickHash> b) {
    const std::size_t shared = std::min(a.size(), b.size());
    for (std::size_t i = 0; i < shared; ++i)
        if (a[i].tick != b[i].tick || a[i].hash != b[i].hash)
            return i;
    if (a.size() != b.size())
        return shared;
    return std::nullopt;
}

double percentile_ms(std::vector<double> values, double percentile);
Result run(Scenario&, const Options&);

} // namespace sqp::benchmark
