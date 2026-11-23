// SPDX-FileCopyrightText: 2026 Erin Catto
// SPDX-License-Identifier: MIT

// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <span>
#include <vector>

#include "sqphysics2d/detail/math.h"
#include "sqphysics2d/settings.h"

namespace sqp::BroadPhase {

// A proxy participates in dynamic or immutable-static pair filtering.
enum ProxyFlags : u8 {
    PROXY_DYNAMIC = 1u << 0,
    PROXY_STATIC = 1u << 1,
    PROXY_KINEMATIC = 1u << 2,
};

struct Proxy {
    u32 id = 0;
    Aabb aabb;
    u32 flags = PROXY_DYNAMIC;
};

// Indices into the broad-phase body/shape columns, not body handles.
struct Pair {
    u32 a = 0;
    u32 b = 0;

    bool operator<(const Pair& rhs) const {
        return a != rhs.a ? a < rhs.a : b < rhs.b;
    }

    bool operator==(const Pair& rhs) const = default;
};

struct QueryStats {
    u64 nodes_visited = 0;
    u64 leaves_visited = 0;
    u64 candidates_emitted = 0;
};

// Cast from origin along translation, clipped to max_fraction of that translation.
struct RayCastInput {
    Vec2 origin;
    Vec2 translation;
    real max_fraction = real(1);
};

struct RayCastStats {
    u64 nodes_visited = 0;
    u64 leaves_visited = 0;
};

// Return a negative value to skip a leaf, zero to terminate, a positive
// fraction less or equal to clipped_input.max_fraction to clip the ray.
using RayCastCallback = real (*)(const RayCastInput& clipped_input, u32 proxy_id, void* context);

} // namespace sqp::BroadPhase

namespace sqp {

namespace detail { class StateHash; }

// Scalar dynamic AABB tree. Proxy handles remain stable across rotations and rebuilds.
class DynamicAabbTree {
public:
    void append_state_hash(detail::StateHash&) const;
    struct QualityStats {
        u32 height = 0;
        u32 reachable_nodes = 0;
        u32 leaf_count = 0;
        double area_ratio = 0;
        double sibling_overlap_ratio = 0;
    };

    std::size_t proxy_count() const { return proxy_count_; }
    // Inspect the three trees in static, kinematic, dynamic order.
    // Returns false if links, heights, or proxy ownership are inconsistent.
    bool inspect(std::array<QualityStats, 3>& quality) const;

    void clear();
    // Rebuilds conservative bounds from exact leaves while preserving proxy handles.
    void set_proxy_extension_ratio(real ratio);
    i32 add_proxy(u32 id, const Aabb& Aabb,
                  u32 flags = static_cast<u32>(BroadPhase::ProxyFlags::PROXY_DYNAMIC));
    i32 add_proxy_unbuilt(u32 id, const Aabb& Aabb,
                          u32 flags = static_cast<u32>(BroadPhase::ProxyFlags::PROXY_DYNAMIC));

    bool update_proxy(i32 handle, const Aabb& Aabb, u32 flags, bool* reinserted = nullptr);
    bool remove_proxy(i32 handle);
    // True when handle is a live leaf registered with the given proxy ID.
    bool proxy_matches(i32 handle, u32 id) const;

    // Rebuild tree links while retaining live proxy handles.
    void rebuild();
    // Box2D-style partial rebuild of changed dynamic and kinematic branches.
    void rebuild_changed();

    void build_pairs(std::vector<BroadPhase::Pair>& out_pairs,
                     std::vector<i32>* traversal_stack = nullptr,
                     u64* ordering_sorts = nullptr) const;

    void query_pairs(u32 id, const Aabb& query, u32 flags, std::vector<BroadPhase::Pair>& out_pairs,
                     std::vector<i32>* traversal_stack = nullptr,
                     BroadPhase::QueryStats* stats = nullptr,
                     std::span<const u32> canonical_query_marks = {},
                     u32 canonical_query_epoch = 0) const;
    // Queries using the bounds, ID and flags already stored for an existing proxy.
    // Returns false for an invalid or removed handle, leaving out_pairs unchanged.
    bool query_proxy_pairs(i32 handle, std::vector<BroadPhase::Pair>& out_pairs,
                           std::vector<i32>* traversal_stack = nullptr,
                           BroadPhase::QueryStats* stats = nullptr,
                           std::span<const u32> canonical_query_marks = {},
                           u32 canonical_query_epoch = 0) const;
    // Append matching proxy IDs to out_ids.
    void query_aabb(const Aabb& query, std::vector<u32>& out_ids) const;
    BroadPhase::RayCastStats ray_cast(const BroadPhase::RayCastInput& input,
                                      BroadPhase::RayCastCallback callback, void* context,
                                      std::vector<i32>* traversal_stack = nullptr) const;

private:
    // Based on Box2D's scalar dynamic tree: siblings occupy adjacent nodes,
    // while stable proxy identities are stored separately from relocatable nodes.
    struct Node {
        Aabb bounds;
        u64 padding = 0;
        u32 flag_index = 0;
        i32 height_or_shape = 0;
    };
    struct Tree {
        std::vector<Node> nodes;
        std::vector<i32> parents;
        std::vector<i32> free_pairs;
        bool ordered = true;
    };
    struct ProxyRecord {
        Aabb exact;
        u32 id = 0;
        u32 flags = 0;
        i32 node = -1;
        u8 tree = 0;
        bool live = false;
    };
    struct BuildItem {
        Aabb bounds;
        i32 node = -1;
        i32 proxy = -1;
    };

    static u8 tree_type(u32 flags);
    static bool is_leaf(const Node& node);
    static bool is_empty(const Node& node);
    static bool moved(const Node& node);
    static i32 child(const Node& node);
    static i32 proxy_index(const Node& node);
    static Node empty_node();
    Node leaf_node(i32 handle, const Aabb& bounds, bool mark) const;
    Node internal_node(const Tree& tree, i32 pair) const;
    void link_children(Tree& tree, i32 node_index);
    void mark_path(Tree& tree, i32 node_index);
    i32 allocate_pair(Tree& tree);
    i32 find_best_sibling(const Tree& tree, const Aabb& bounds) const;
    void rotate(Tree& tree, i32 index);
    void insert_leaf(Tree& tree, i32 handle, const Aabb& bounds, bool rotate_nodes);
    void remove_leaf(Tree& tree, i32 handle);
    void gather_build_items(const Tree& tree, i32 node, std::vector<BuildItem>& items,
                            bool full) const;
    void copy_subtree(const Tree& old_tree, i32 old_node, Tree& new_tree, i32 new_node,
                      i32 parent);
    void build_subtree(const Tree& old_tree, std::vector<BuildItem>& items,
                       std::size_t first, std::size_t last, Tree& new_tree, i32 target,
                       i32 parent);
    void rebuild_tree(u8 type, bool full);
    bool valid_handle(i32 handle) const;
    bool should_pair(u32 a, u32 b) const;
    void query_tree(u8 type, u32 id, const Aabb& bounds, u32 flags,
                    std::vector<BroadPhase::Pair>& out, std::vector<i32>* scratch,
                    BroadPhase::QueryStats* stats,
                    std::span<const u32> marks, u32 epoch) const;

    real proxy_extension_ratio_ = BroadPhaseSettings::default_proxy_extension_ratio;
    std::array<Tree, 3> trees_;
    // Reused only while rebuilding; the active trees retain their original layout.
    Tree rebuild_scratch_;
    std::vector<BuildItem> rebuild_items_;
    std::vector<ProxyRecord> proxies_;
    std::size_t proxy_count_ = 0;
};

} // namespace sqp
