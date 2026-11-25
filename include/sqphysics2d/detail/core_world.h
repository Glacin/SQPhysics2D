// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#pragma once

#include "broad_phase.h"
#include "constraint_graph.h"
#include "contact.h"
#include "contact_solver.h"
#include "data_storage.h"
#include "shape_query_view.h"
#include "sqphysics2d/events.h"
#include "sqphysics2d/joint.h"
#include "sqphysics2d/settings.h"
#include "sqphysics2d/world.h"

#include <array>
#include <cassert>
#include <cstdio>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace sqp {

// Dispatches stage functions until one returns Finish or Error.
// Each stage chooses the next stage ID; this base class holds no physics policy.
class VirtualWorld {
public:
    using StageId = u32;
    using StageFn = StageId (*)(VirtualWorld&);

    static constexpr StageId Finish = 0;
    static constexpr StageId Error = 1;

    virtual ~VirtualWorld() = default;

    // Executes stages until Finish or Error.
    bool step() {
        last_error_message_.clear();
        StageId id = entry_stage_;
        StageId previous = Finish;

        while (id != Finish) {
            if (id == Error || id >= stages_.size() || !stages_[id]) {
                const std::string source = previous != Finish && previous < stage_names_.size()
                                               ? stage_names_[previous]
                                               : "dispatcher";
                last_error_message_ = id == Error
                                          ? source + ": stage returned Error"
                                          : source + ": invalid stage ID " + std::to_string(id);
                std::fprintf(stderr, "%s\n", last_error_message_.c_str());
                return false;
            }
            previous = id;
            id = stages_[id](*this);
        }
        return true;
    }

    std::string_view last_error_message() const noexcept { return last_error_message_; }

protected:
    // Registers built-in terminal stages.
    VirtualWorld() {
        stages_ = {nullptr, nullptr};
        stage_names_ = {"Finish", "Error"};
    }

    // Registers a stage and returns its assigned ID.
    StageId install_stage(StageFn fn, std::string name) {
        assert(fn);
        assert(stages_.size() < std::numeric_limits<StageId>::max());

        StageId id = static_cast<StageId>(stages_.size());
        stage_names_.push_back(std::move(name));
        stages_.push_back(fn);
        return id;
    }

    void set_entry_stage(StageId id) {
        entry_stage_ = id;
    }

    std::vector<StageFn> stages_;
    std::vector<std::string> stage_names_;
    StageId entry_stage_ = Finish;
    std::string last_error_message_;
};

// Simulation ownership and scratch storage for the default stage pipeline.
// Stage vectors are reused between intervals; consumers must observe each stage's alignment rules.
struct WorldState {
    WorldState() = default;
    WorldState(const WorldState&);
    WorldState& operator=(const WorldState&) = delete;
    WorldState(WorldState&&) = delete;
    WorldState& operator=(WorldState&&) = delete;
    // The first failing check supplies detail; the dispatcher adds the stage name.
    std::string stage_error_detail;
    std::string last_error_message;
    struct StageStats {
        enum class Stage : u8 {
            ProcessCommand,
            CollectTasks,
            BroadPhase,
            NarrowPhase,
            ProcessToi,
            UpdatePairs,
            PreSolve,
            Solve,
            PostSolve,
            AdvanceInterval,
            Count,
        };
        static constexpr std::size_t stage_count = static_cast<std::size_t>(Stage::Count);
        struct Entry {
            u64 calls = 0;
            u64 nanoseconds = 0;
        };
        struct IntervalEntry {
            u64 active_islands = 0;
            u64 active_bodies = 0;
            u64 candidates = 0;
            u64 pairs = 0;
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
            u64 query_nanoseconds = 0;
            u64 update_nanoseconds = 0;
            u64 rebuild_nanoseconds = 0;
            u64 nodes_visited = 0;
            u64 proxy_reinsertions = 0;
            u64 query_calls = 0;
#endif
        };
        u64 intervals = 0;
        struct BodyFreeze {
            BodyHandle body;
            tick_t start = 0, end = 0;
            bool exhausted = false;
        };
        bool record_body_freezes = false;
        std::vector<BodyFreeze> body_freezes;
        struct ToiOutcomes {
            u64 microsteps = 0, effective_microsteps = 0;
            u64 islands = 0, charged_islands = 0;
            u64 selected_pairs = 0, solved_pairs = 0, effective_pairs = 0;
            u64 approaching_contact_pairs = 0;
        } toi_outcomes;
        bool record_toi_outcomes = false;
        std::array<Entry, stage_count> stages{};
#ifdef SQP_ENABLE_PERFORMANCE_PROFILING
        u64 pending_candidates = 0;
        u64 pending_query_nanoseconds = 0;
        u64 pending_update_nanoseconds = 0;
        u64 pending_rebuild_nanoseconds = 0;
        u64 pending_nodes_visited = 0;
        u64 pending_proxy_reinsertions = 0;
        u64 pending_query_calls = 0;
        u64 freeze_retries = 0;
        u64 freeze_retry_queries = 0;
        u64 freeze_retry_sweeps = 0;
        u64 narrow_sweeps = 0;
        u64 narrow_dimension_requests = 0;
        u64 narrow_dimension_loads = 0;
        u64 narrow_start_predictions = 0;
        u64 narrow_reused_start_predictions = 0;
        u64 narrow_reused_results = 0;
        u64 toi_dedup_checks = 0;
        u64 toi_duplicates = 0;
        u64 constraint_prepare_nanoseconds = 0;
        std::vector<IntervalEntry> interval_entries;
#endif
    } stats;

    PhysicsSettings settings;
    BodyPool body_pool;
    ShapePool shape_pool;
    IslandPool island_pool;
    JointPool joint_pool;
    ConstraintGraph contact_graph;

    struct ConnectivityData {
        u64 seen_graph_version = std::numeric_limits<u64>::max();
        bool dirty = true;
        bool full_rebuild = true;
        std::vector<IslandId> dirty_islands;
        std::size_t last_rebuilt_body_count = 0;
        std::vector<BodyHandle> bodies;
        std::vector<u32> parents;
        std::vector<u32> roots;
        std::vector<IslandId> old_islands;
        std::vector<IslandId> destinations;
        std::vector<u32> component_toi_truncations;
        std::vector<u8> claimed_slots;
        std::array<std::vector<u32>, (size_t)(BodyConfigs::Kind::AabbBox) + 1> body_positions;
        std::array<std::vector<u32>, (size_t)(BodyConfigs::Kind::AabbBox) + 1> position_stamps;
        u32 position_epoch = 0;
    } connectivity;

    struct CommandData {
        // The heap orders absolute times, then command IDs; results are consumed by ID.
        struct Pending {
            tick_t time = 0;
            CommandId id = 0;
            PhysicsCommand command;
        };
        struct Later {
            bool operator()(const Pending& a, const Pending& b) const {
                return a.time != b.time ? a.time > b.time : a.id > b.id;
            }
        };
        struct Completed {
            CommandId id = 0;
            CommandResult result;
        };
        std::vector<Pending> heap;
        std::vector<Completed> results;
        CommandId next_id = 1;
    } commands;

    struct BroadPhaseData {
        // Interval data of the broad-phase stage.
        DynamicAabbTree tree;
        std::vector<BodyHandle> bodies;
        std::array<std::vector<std::vector<u32>>, (size_t)(BodyConfigs::Kind::AabbBox) + 1>
            proxies_by_dense;
        bool proxy_index_ready = false;
        std::vector<ShapeHandle> shapes;
        std::vector<Aabb> start_aabbs;
        std::vector<Aabb> end_aabbs;
        std::vector<i32> proxy_handles;
        // A missing proxy is valid only when this entry was suppressed by NoCollision.
        std::vector<u8> proxy_suppressed;

        // Reused configuration and candidate buffers for stage-level filtering.
        std::vector<BodyConfiguration> filter_configurations;
        std::vector<BroadPhase::Pair> candidate_query_scratch;

        // Stack for tree traversal during pair generation and queries.
        std::vector<i32> traversal_stack;

        // Indices of proxies updated in this interval; reused for pair queries.
        std::vector<u32> interval_proxy_indices;
        std::vector<u32> selected;
        u32 selection_epoch = 0;
        // Islands whose proxy sweeps changed during TOI solving or freeze propagation.
        std::vector<IslandId> toi_changed_islands;
        bool reprepare_after_toi = false;
        bool freeze_retry = false;

        // Output of BroadPhase stage.
        std::vector<BroadPhase::Pair> candidates;
        // Candidate set for the ordinary interval, reused after a TOI microstep.
        std::vector<BroadPhase::Pair> ordinary_candidates;
    } broad_phase;

    // Raw sweep impact pose in the pair's frame, anchored at A's sweep start.
    struct ToiLocalImpact {
        Vec2 origin;
        Vec2 position_a, position_b;
        real rotation_a = 0, rotation_b = 0;
        bool valid = false;
    };

    // Solver geometry retained before world-space manifold point rounding.
    struct ToiLocalContact {
        Vec2 origin;
        ContactManifold manifold;
        bool valid = false;
    };

    struct NarrowPhaseData {

        struct Pose {
            Vec2 position;
            real rotation = 0;
        };

        struct CandidateCacheEntry {
            u32 epoch = 0;
            bool valid = false;
            BodyHandle body;
            ShapeHandle shape;
            detail::ShapeQueryView view;
            Pose source, end, predicted_start;
            tick_t predicted_start_time = 0;
            bool frozen = false;
            real sweep_radius = 0;
            bool dimensions_loaded = false;
            real radius = 0;
            Vec2 half_extents;
            Vec2 linear_velocity;
            real angular_velocity = 0;
            BodyConfiguration configuration;
            tick_t island_time = 0;
            u64 motion_version = 0;
        };

        // Indexed by broad_phase.bodies; entries are valid only for the current stage run.
        std::vector<CandidateCacheEntry> candidate_cache;
        u32 candidate_cache_epoch = 0;

        struct SweepHit {
            tick_t time;
            Vec2 normal;
            ToiLocalImpact impact;
        };
        struct PairResult {
            BroadPhase::Pair candidate;
            Pair pair;
            tick_t start_time = 0;
            u8 modes = 0;
            u8 frozen = 0;
            // Absolute endpoints, followed by start poses and velocities for a sweep.
            std::array<real, 18> motion{};
            std::variant<std::monostate, ContactManifold, SweepHit> result;
        };
        std::vector<PairResult> pair_results, next_pair_results;
        tick_t pair_result_end = 0;

        // Ordinary contact rows are reusable only within the same prepared time span.
        std::optional<tick_t> prepared_time;
        tick_t prepared_end = 0;
        std::vector<BroadPhase::Pair> ordinary_contact_keys;
        std::vector<BroadPhase::Pair> previous_contact_keys;

        // Aligned TOI event columns; toi_impacts records the unrounded impact pose.
        std::vector<Pair> toi_pairs;
        // Absolute TOI time.
        std::vector<tick_t> toi_times;
        // TOI normal directed from B toward A.
        std::vector<Vec2> toi_normals;
        std::vector<ToiLocalImpact> toi_impacts;

    } narrow_phase;

    struct ProcessToiData {

        struct ToiEvent {
            tick_t time = 0;
            Pair pair;
            Vec2 normal;
            u64 sequence = 0;
            ToiLocalImpact impact;
        };
        struct ToiEventCompare {
            bool operator()(const ToiEvent& a, const ToiEvent& b) const {
                if (a.time != b.time)
                    return a.time > b.time;
                return a.sequence > b.sequence;
            }
        };

        struct EventKey {
            tick_t time;
            Pair pair;
            bool operator==(const EventKey&) const = default;
        };
        struct EventKeyHash {
            std::size_t operator()(const EventKey& key) const noexcept {
                using Hash = detail::ContactPairHash;
                const u64 shapes = Hash::mix(key.pair.shape_a.index, key.pair.shape_b.index);
                return static_cast<std::size_t>(Hash::mix(shapes, key.time));
            }
        };

        // Future events survive interval changes and are revalidated before use.
        std::vector<ToiEvent> toi_heap;
        std::unordered_set<EventKey, EventKeyHash> toi_keys;
        u64 next_sequence = 0;
        std::vector<ToiEvent> selected_events;

    } process_toi;

    struct PairsData {
        // Active pair, world manifold, and optional pair-local manifold rows.
        std::unique_ptr<std::vector<Pair>> pairs = std::make_unique<std::vector<Pair>>();
        std::unique_ptr<ContactManifolds> manifolds = std::make_unique<ContactManifolds>();
        std::vector<ToiLocalContact> toi_local_contacts;
        // NarrowPhase input retained across PreSolve changes.
        ContactManifolds original_manifolds;
        bool original_manifolds_valid = false;

        // Ordinary results parked while the active columns hold microstep results.
        std::unique_ptr<std::vector<Pair>> delayed_pairs = std::make_unique<std::vector<Pair>>();
        std::unique_ptr<ContactManifolds> delayed_manifolds = std::make_unique<ContactManifolds>();
        std::vector<ToiLocalContact> delayed_toi_local_contacts;
        ContactManifolds delayed_original_manifolds;
        bool delayed_original_manifolds_valid = false;

        struct PreSolveScratch {
            ContactManifolds manifolds;
            std::vector<ToiLocalContact> local_contacts;
            std::vector<real> restitution, dynamic_friction, static_friction;
            std::vector<real> damping, frequency, correction_speed;
        } pre_solve_scratch;

        // PreSolve material output; every column has the same row index as pairs.
        std::vector<real> combined_damping_ratios;
        std::vector<real> combined_contact_frequencies;
        std::vector<real> combined_max_penetration_correction_speeds;
        std::vector<real> combined_static_frictions;
        std::vector<real> combined_dynamic_frictions;
        std::vector<real> combined_restitutions;

    } pairs;

    struct IntervalData {
        // Active islands for the current ordinary or microstep interval.
        std::unique_ptr<std::vector<IslandId>> islands = std::make_unique<std::vector<IslandId>>();
        // Ordinary islands parked while the active list holds microstep islands.
        std::unique_ptr<std::vector<IslandId>> delayed_islands =
            std::make_unique<std::vector<IslandId>>();

        bool microstep = false;
        bool tick = false;
        tick_t microstep_time = 0;

        tick_t target_end = 0;
    } interval;

    struct ScheduleData {
        struct Entry {
            IslandId id;
            tick_t period = 0;
            tick_t due = 0;
            std::size_t position = 0;
            // Join earlier ordinary intervals chosen for other islands to keep CCD sweeps current.
            bool uses_ccd = false;
        };
        // One bucket per power-of-two period; entries are indexed by Island slot.
        // Bucket zero is parked sleep; active periods use countr_zero(period) + 1.
        std::array<std::vector<IslandId>, config::PERIOD_BUCKET_COUNT + 1> period_islands;
        std::vector<Entry> entries;
        std::vector<IslandId> scratch_islands;
        tick_t time = 0;
        tick_t step_target = 0;
        tick_t ordinary_target = 0;
        std::optional<tick_t> last_ordinary_end;
        bool prepared = false;
    } schedule;

    // Reusable solve selection, constraint and body-integration buffers.
    struct SolveData {
        struct ContactCacheEntry {
            BodyHandle body_a, body_b;
            ShapeHandle shape_a, shape_b;
            ContactPointKey keys[2];
            real normal_impulses[2]{}, tangent_impulses[2]{};
            Vec2 normal;
            real step_seconds = 0;
            tick_t time = 0;
            u8 count = 0;
        };
        // Persists independently of the interval's active and delayed pair columns.
        std::vector<ContactCacheEntry> contact_cache;
        std::unordered_map<detail::ContactPairKey, std::size_t, detail::ContactPairHash>
            contact_cache_index;
        std::vector<std::size_t> constraint_cache_slots;
        struct ProcessedIsland {
            IslandId id;
            tick_t start = 0;
            tick_t end = 0;
            bool received_toi_impulse = false;
        };
        std::vector<ProcessedIsland> processed_islands;
        // Rebuilt each invocation; storage capacity is retained across intervals.
        std::vector<IslandId> islands;
        std::vector<BodyHandle> bodies;
        std::vector<tick_t> body_start_times;

        // Contact-connected islands share a solve start. Indexed by selected island.
        std::vector<std::size_t> component_parents;
        std::vector<std::size_t> island_positions;
        std::vector<u32> island_position_stamps;
        u32 island_position_epoch = 0;
        std::vector<tick_t> component_starts;
        std::vector<real> body_step_seconds;
        std::vector<real> constraint_step_seconds;

        // Indexed by island slot and by body kind/dense index, respectively.
        // An epoch avoids clearing these buffers on every Solve invocation.
        std::vector<u32> island_stamps;
        std::array<std::vector<u32>, static_cast<std::size_t>(BodyConfigs::Kind::AabbBox) + 1>
            body_stamps;
        u32 selection_epoch = 0;
        std::vector<ContactConstraint> constraints;
        std::vector<PreparedContactConstraint> prepared;
        std::vector<std::size_t> constraint_pair_indices;

    } solve;

    struct PostSolveData {
        struct Report {
            std::size_t pair_index = 0;
            ContactSolveReport contact;
        };
        // Only solved pairs with at least one AfterSolve subscriber are retained.
        std::vector<Report> reports;
        std::vector<PostSolveEvent> events;
    } post_solve;

    struct CallbackData {
        // Registered callbacks and their opaque user pointers.
        EventCallback pre_solve = nullptr;
        EventCallback post_solve = nullptr;
        void* pre_data = nullptr;
        void* post_data = nullptr;
    } callbacks;
};

// Installs the default physics stages and owns the state accessed by World.
// Stage entry points return the ID of the next stage to dispatch.
class CoreWorld : public VirtualWorld {
public:
    ~CoreWorld() override;

protected:
    friend class World;
    CoreWorld();
    explicit CoreWorld(const WorldState& source);
    void install_default_stages();

    bool step_until(tick_t target);
    bool process_commands();

    enum class StageRole : u8 {
        ProcessCommand,
        CollectTasks,
        BroadPhase,
        NarrowPhase,
        ProcessToi,
        UpdatePairs,
        PreSolve,
        Solve,
        PostSolve,
        // Choose next interval target time.
        AdvanceInterval,
        Count,
    };

    // Default physics stage implementations.
    static StageId process_command(VirtualWorld& world);
    static StageId collect_tasks(VirtualWorld& world);
    static StageId broad_phase(VirtualWorld& world);
    static StageId update_pairs(VirtualWorld& world);
    static StageId narrow_phase(VirtualWorld& world);
    static StageId process_toi(VirtualWorld& world);
    static StageId pre_solve(VirtualWorld& world);
    static StageId solve(VirtualWorld& world);
    static StageId post_solve(VirtualWorld& world);
    static StageId advance_interval(VirtualWorld& world);

    StageId stage_id(StageRole role) const {
        return stage_ids_[static_cast<size_t>(role)];
    }

    std::array<StageId, static_cast<size_t>(StageRole::Count)> stage_ids_{};

    // Physics state and working data.
    WorldState world_state_;
};

} // namespace sqp
