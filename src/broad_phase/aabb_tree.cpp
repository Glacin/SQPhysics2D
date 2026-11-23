// SPDX-FileCopyrightText: 2026 Erin Catto
// SPDX-License-Identifier: MIT
//
// Scalar sibling-pair layout, insertion, movement and partial rebuilding adapted
// from the vendored Box2D dynamic_tree.c (Box2D 3.2.0).
// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include "sqphysics2d/detail/broad_phase.h"
#include "world/state_hash.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace sqp {
namespace {
constexpr u32 moved_bit = 1u << 30;
constexpr u32 leaf_bit = 1u << 31;
constexpr u32 index_mask = ~(moved_bit | leaf_bit);
constexpr u32 empty_bits = index_mask | leaf_bit;

real perimeter(const Aabb& a) {
    return real(2) * (a.max.x - a.min.x + a.max.y - a.min.y);
}
Aabb expanded_bounds(const Aabb& exact, u32 flags, real ratio) {
    if (flags & BroadPhase::PROXY_STATIC) return exact;
    const Vec2 e = exact.extent() * ratio;
    return {exact.min - e, exact.max + e};
}
bool contains(const Aabb& outer, const Aabb& inner) {
    return outer.min.x <= inner.min.x && outer.min.y <= inner.min.y &&
           outer.max.x >= inner.max.x && outer.max.y >= inner.max.y;
}
real distance_squared(Vec2 a, Vec2 b) {
    const Vec2 d = a - b;
    return math::dot(d, d);
}
} // namespace

u8 DynamicAabbTree::tree_type(u32 flags) {
    return flags & BroadPhase::PROXY_STATIC ? 0 :
           flags & BroadPhase::PROXY_KINEMATIC ? 1 : 2;
}
bool DynamicAabbTree::is_leaf(const Node& n) { return (n.flag_index & leaf_bit) != 0; }
bool DynamicAabbTree::is_empty(const Node& n) { return n.flag_index == empty_bits; }
bool DynamicAabbTree::moved(const Node& n) { return (n.flag_index & moved_bit) != 0; }
i32 DynamicAabbTree::child(const Node& n) { return static_cast<i32>(n.flag_index & index_mask); }
i32 DynamicAabbTree::proxy_index(const Node& n) {
    return static_cast<i32>(n.flag_index & index_mask);
}
DynamicAabbTree::Node DynamicAabbTree::empty_node() {
    Node n;
    n.flag_index = empty_bits;
    return n;
}
DynamicAabbTree::Node DynamicAabbTree::leaf_node(i32 handle, const Aabb& bounds,
                                                  bool mark) const {
    Node n;
    n.bounds = bounds;
    n.flag_index = static_cast<u32>(handle) | leaf_bit | (mark ? moved_bit : 0);
    n.height_or_shape = static_cast<i32>(proxies_[handle].id);
    return n;
}
DynamicAabbTree::Node DynamicAabbTree::internal_node(const Tree& tree, i32 pair) const {
    const Node& a = tree.nodes[pair];
    const Node& b = tree.nodes[pair + 1];
    Node n;
    n.bounds = Aabb::merge(a.bounds, b.bounds);
    n.flag_index = static_cast<u32>(pair) |
                   ((a.flag_index | b.flag_index) & moved_bit);
    n.height_or_shape = 1 + std::max(is_leaf(a) ? 0 : a.height_or_shape,
                                     is_leaf(b) ? 0 : b.height_or_shape);
    return n;
}
void DynamicAabbTree::link_children(Tree& tree, i32 index) {
    const Node& n = tree.nodes[index];
    if (is_leaf(n)) {
        proxies_[proxy_index(n)].node = index;
    } else {
        const i32 pair = child(n);
        tree.parents[pair] = index;
        tree.parents[pair + 1] = index;
    }
}
void DynamicAabbTree::mark_path(Tree& tree, i32 index) {
    while (index >= 0) {
        tree.nodes[index].flag_index |= moved_bit;
        index = tree.parents[index];
    }
}
i32 DynamicAabbTree::allocate_pair(Tree& tree) {
    if (!tree.free_pairs.empty()) {
        const i32 pair = tree.free_pairs.back();
        tree.free_pairs.pop_back();
        return pair;
    }
    const i32 pair = static_cast<i32>(tree.nodes.size());
    tree.nodes.resize(tree.nodes.size() + 2);
    tree.parents.resize(tree.parents.size() + 2, -1);
    return pair;
}
i32 DynamicAabbTree::find_best_sibling(const Tree& tree, const Aabb& leaf) const {
    const auto& nodes = tree.nodes;
    i32 index = 0, best = 0;
    const Vec2 center = leaf.center();
    const real area = perimeter(leaf);
    real area_base = perimeter(nodes[0].bounds);
    real direct = perimeter(Aabb::merge(nodes[0].bounds, leaf));
    real inherited = 0, best_cost = direct;
    while (!is_leaf(nodes[index])) {
        const real cost = direct + inherited;
        if (cost < best_cost) { best = index; best_cost = cost; }
        inherited += direct - area_base;
        const i32 pair = child(nodes[index]);
        real lower[2] = {std::numeric_limits<real>::max(),
                         std::numeric_limits<real>::max()};
        real child_area[2] = {}, child_direct[2] = {};
        bool leaf_child[2] = {};
        for (i32 side = 0; side < 2; ++side) {
            const i32 c = pair + side;
            leaf_child[side] = is_leaf(nodes[c]);
            child_direct[side] = perimeter(Aabb::merge(nodes[c].bounds, leaf));
            if (leaf_child[side]) {
                const real candidate = child_direct[side] + inherited;
                if (candidate < best_cost) { best = c; best_cost = candidate; }
            } else {
                child_area[side] = perimeter(nodes[c].bounds);
                lower[side] = inherited + child_direct[side] +
                              std::min(area - child_area[side], real(0));
            }
        }
        if ((leaf_child[0] && leaf_child[1]) ||
            (best_cost <= lower[0] && best_cost <= lower[1])) break;
        if (lower[0] == lower[1] && !leaf_child[0] && !leaf_child[1]) {
            lower[0] = distance_squared(nodes[pair].bounds.center(), center);
            lower[1] = distance_squared(nodes[pair + 1].bounds.center(), center);
        }
        const i32 side = lower[0] < lower[1] && !leaf_child[0] ? 0 : 1;
        if (leaf_child[side]) break;
        index = pair + side;
        area_base = child_area[side];
        direct = child_direct[side];
    }
    return best;
}
void DynamicAabbTree::rotate(Tree& tree, i32 index) {
    const i32 pair = child(tree.nodes[index]);
    const i32 left = pair, right = pair + 1;
    if (is_leaf(tree.nodes[left]) && is_leaf(tree.nodes[right])) return;
    i32 down = -1, up = -1;
    real best_delta = 0;
    if (!is_leaf(tree.nodes[right])) {
        const i32 c = child(tree.nodes[right]);
        const real old = perimeter(tree.nodes[right].bounds);
        const real d0 = perimeter(Aabb::merge(tree.nodes[left].bounds,
                                               tree.nodes[c + 1].bounds)) - old;
        const real d1 = perimeter(Aabb::merge(tree.nodes[left].bounds,
                                               tree.nodes[c].bounds)) - old;
        if (d0 < best_delta) { down = left; up = c; best_delta = d0; }
        if (d1 < best_delta) { down = left; up = c + 1; best_delta = d1; }
    }
    if (!is_leaf(tree.nodes[left])) {
        const i32 c = child(tree.nodes[left]);
        const real old = perimeter(tree.nodes[left].bounds);
        const real d0 = perimeter(Aabb::merge(tree.nodes[right].bounds,
                                               tree.nodes[c + 1].bounds)) - old;
        const real d1 = perimeter(Aabb::merge(tree.nodes[right].bounds,
                                               tree.nodes[c].bounds)) - old;
        if (d0 < best_delta) { down = right; up = c; best_delta = d0; }
        if (d1 < best_delta) { down = right; up = c + 1; best_delta = d1; }
    }
    if (down < 0) return;
    std::swap(tree.nodes[down], tree.nodes[up]);
    link_children(tree, down);
    link_children(tree, up);
    const i32 sibling = down ^ 1;
    tree.nodes[sibling] = internal_node(tree, child(tree.nodes[sibling]));
    tree.ordered = false;
}
void DynamicAabbTree::insert_leaf(Tree& tree, i32 handle, const Aabb& bounds,
                                   bool rotate_nodes) {
    if (tree.nodes.empty()) {
        tree.nodes.resize(2, empty_node());
        tree.parents.resize(2, -1);
    }
    const Node leaf = leaf_node(handle, bounds, true);
    if (is_empty(tree.nodes[0])) {
        tree.nodes[0] = leaf;
        proxies_[handle].node = 0;
        return;
    }
    const i32 sibling = find_best_sibling(tree, bounds);
    const i32 pair = allocate_pair(tree);
    tree.nodes[pair] = tree.nodes[sibling];
    tree.nodes[pair + 1] = leaf;
    tree.parents[pair] = tree.parents[pair + 1] = sibling;
    link_children(tree, pair);
    proxies_[handle].node = pair + 1;
    tree.nodes[sibling] = internal_node(tree, pair);
    tree.ordered = false;
    for (i32 index = sibling; index >= 0; index = tree.parents[index]) {
        if (rotate_nodes) rotate(tree, index);
        tree.nodes[index] = internal_node(tree, child(tree.nodes[index]));
    }
}
void DynamicAabbTree::remove_leaf(Tree& tree, i32 handle) {
    const i32 leaf = proxies_[handle].node;
    if (leaf == 0) {
        tree.nodes[0] = empty_node();
        proxies_[handle].node = -1;
        return;
    }
    const i32 parent = tree.parents[leaf];
    const i32 ancestor = tree.parents[parent];
    tree.nodes[parent] = tree.nodes[leaf ^ 1];
    link_children(tree, parent);
    tree.free_pairs.push_back(leaf & ~1);
    proxies_[handle].node = -1;
    mark_path(tree, parent);
    for (i32 index = ancestor; index >= 0; index = tree.parents[index])
        tree.nodes[index] = internal_node(tree, child(tree.nodes[index]));
}
void DynamicAabbTree::clear() {
    trees_ = {};
    proxies_.clear();
    proxy_count_ = 0;
}
bool DynamicAabbTree::valid_handle(i32 handle) const {
    return handle >= 0 && static_cast<std::size_t>(handle) < proxies_.size() &&
           proxies_[handle].live;
}
void DynamicAabbTree::set_proxy_extension_ratio(real ratio) {
    ratio = std::isfinite(ratio) ? std::clamp(ratio, real(0), real(1))
                                 : BroadPhaseSettings::default_proxy_extension_ratio;
    if (ratio == proxy_extension_ratio_) return;
    proxy_extension_ratio_ = ratio;
    rebuild();
}
i32 DynamicAabbTree::add_proxy(u32 id, const Aabb& bounds, u32 flags) {
    if (!bounds.is_valid() || proxies_.size() >= index_mask) return -1;
    const i32 handle = static_cast<i32>(proxies_.size());
    const u8 type = tree_type(flags);
    proxies_.push_back({bounds, id, flags, -1, type, true});
    ++proxy_count_;
    insert_leaf(trees_[type], handle, expanded_bounds(bounds, flags, proxy_extension_ratio_),
                true);
    return handle;
}
i32 DynamicAabbTree::add_proxy_unbuilt(u32 id, const Aabb& bounds, u32 flags) {
    if (!bounds.is_valid() || proxies_.size() >= index_mask) return -1;
    const i32 handle = static_cast<i32>(proxies_.size());
    proxies_.push_back({bounds, id, flags, -1, tree_type(flags), true});
    ++proxy_count_;
    return handle;
}
bool DynamicAabbTree::update_proxy(i32 handle, const Aabb& bounds, u32 flags, bool* reinserted) {
    if (reinserted) *reinserted = false;
    if (!valid_handle(handle) || !bounds.is_valid()) return false;
    auto& record = proxies_[handle];
    const u8 type = tree_type(flags);
    const bool type_changed = type != record.tree;
    if (record.exact == bounds && !type_changed) {
        record.flags = flags;
        return true;
    }
    if (record.node >= 0 && !type_changed && type != 0 &&
        contains(trees_[type].nodes[record.node].bounds, bounds)) {
        record.exact = bounds;
        record.flags = flags;
        return true;
    }
    if (record.node >= 0) remove_leaf(trees_[record.tree], handle);
    record.exact = bounds;
    record.flags = flags;
    record.tree = type;
    insert_leaf(trees_[type], handle, expanded_bounds(bounds, flags, proxy_extension_ratio_),
                false);
    if (reinserted) *reinserted = true;
    return true;
}
bool DynamicAabbTree::remove_proxy(i32 handle) {
    if (!valid_handle(handle)) return false;
    auto& record = proxies_[handle];
    if (record.node >= 0) remove_leaf(trees_[record.tree], handle);
    record.live = false;
    --proxy_count_;
    return true;
}
bool DynamicAabbTree::proxy_matches(i32 handle, u32 id) const {
    return valid_handle(handle) && proxies_[handle].id == id;
}
bool DynamicAabbTree::should_pair(u32 a, u32 b) const {
    return !((a & BroadPhase::PROXY_STATIC) && (b & BroadPhase::PROXY_STATIC));
}
void DynamicAabbTree::gather_build_items(const Tree& tree, i32 node,
                                         std::vector<BuildItem>& items, bool full) const {
    const Node& n = tree.nodes[node];
    if (!full && !moved(n)) {
        items.push_back({n.bounds, node, -1});
    } else if (is_leaf(n)) {
        items.push_back({n.bounds, -1, proxy_index(n)});
    } else {
        const i32 pair = child(n);
        gather_build_items(tree, pair, items, full);
        gather_build_items(tree, pair + 1, items, full);
    }
}
void DynamicAabbTree::copy_subtree(const Tree& old_tree, i32 old_node, Tree& new_tree,
                                    i32 new_node, i32 parent) {
    const Node old = old_tree.nodes[old_node];
    new_tree.nodes[new_node] = old;
    new_tree.nodes[new_node].flag_index &= ~moved_bit;
    new_tree.parents[new_node] = parent;
    if (is_leaf(old)) {
        proxies_[proxy_index(old)].node = new_node;
        return;
    }
    const i32 old_pair = child(old);
    const i32 pair = allocate_pair(new_tree);
    copy_subtree(old_tree, old_pair, new_tree, pair, new_node);
    copy_subtree(old_tree, old_pair + 1, new_tree, pair + 1, new_node);
    new_tree.nodes[new_node] = internal_node(new_tree, pair);
}
void DynamicAabbTree::build_subtree(const Tree& old_tree, std::vector<BuildItem>& items,
                                    std::size_t first, std::size_t last, Tree& new_tree,
                                    i32 target, i32 parent) {
    if (last - first == 1) {
        const BuildItem item = items[first];
        if (item.proxy >= 0) {
            new_tree.nodes[target] = leaf_node(item.proxy, item.bounds, false);
            new_tree.parents[target] = parent;
            proxies_[item.proxy].node = target;
        } else {
            copy_subtree(old_tree, item.node, new_tree, target, parent);
        }
        return;
    }
    Aabb centers{items[first].bounds.center(), items[first].bounds.center()};
    for (std::size_t i = first + 1; i < last; ++i)
        centers.include(items[i].bounds.center());
    const bool x = centers.extent().x >= centers.extent().y;
    std::size_t middle = first + (last - first) / 2;
    bool partitioned = false;
    // Box2D's eight-bin SAH rebuild heuristic, evaluated on both axes.
    struct Bin {
        Aabb bounds;
        std::size_t count = 0;
    };
    const auto combine = [](const Bin& a, const Bin& b) {
        if (a.count == 0) return b;
        if (b.count == 0) return a;
        return Bin{Aabb::merge(a.bounds, b.bounds), a.count + b.count};
    };
    real best_cost = std::numeric_limits<real>::max();
    i32 best_axis = -1, best_split = 0;
    real best_min = 0, best_spread = 0;
    for (i32 axis = 0; axis < 2; ++axis) {
        const real minimum = axis == 0 ? centers.min.x : centers.min.y;
        const real spread = axis == 0 ? centers.max.x - minimum : centers.max.y - minimum;
        if (spread <= 0) continue;
        const auto bin_of = [axis, minimum, spread](const BuildItem& item) {
            const real center = axis == 0 ? item.bounds.center().x : item.bounds.center().y;
            return std::clamp(static_cast<i32>(real(8) * (center - minimum) / spread), 0, 7);
        };
        std::array<Bin, 8> bins{};
        for (std::size_t i = first; i < last; ++i) {
            Bin& bin = bins[bin_of(items[i])];
            bin.bounds = bin.count == 0 ? items[i].bounds
                                        : Aabb::merge(bin.bounds, items[i].bounds);
            ++bin.count;
        }
        std::array<Bin, 8> prefix{}, suffix{};
        prefix[0] = bins[0];
        suffix[7] = bins[7];
        for (i32 bin = 1; bin < 8; ++bin)
            prefix[bin] = combine(prefix[bin - 1], bins[bin]);
        for (i32 bin = 6; bin >= 0; --bin)
            suffix[bin] = combine(bins[bin], suffix[bin + 1]);
        for (i32 split = 1; split < 8; ++split) {
            const Bin& left = prefix[split - 1];
            const Bin& right = suffix[split];
            if (left.count == 0 || right.count == 0) continue;
            const real cost = static_cast<real>(left.count) * perimeter(left.bounds) +
                              static_cast<real>(right.count) * perimeter(right.bounds);
            if (cost < best_cost) {
                best_cost = cost;
                best_axis = axis;
                best_split = split;
                best_min = minimum;
                best_spread = spread;
            }
        }
    }
    if (best_axis >= 0) {
        const auto cut = std::partition(items.begin() + first, items.begin() + last,
                                        [best_axis, best_split, best_min, best_spread](
                                            const BuildItem& item) {
            const real center = best_axis == 0 ? item.bounds.center().x
                                               : item.bounds.center().y;
            const i32 bin = std::clamp(
                static_cast<i32>(real(8) * (center - best_min) / best_spread), 0, 7);
            return bin < best_split;
        });
        middle = static_cast<std::size_t>(cut - items.begin());
        partitioned = first < middle && middle < last;
    }
    if (!partitioned)
        std::nth_element(items.begin() + first, items.begin() + middle, items.begin() + last,
                         [x](const BuildItem& a, const BuildItem& b) {
                             const real ca = x ? a.bounds.center().x : a.bounds.center().y;
                             const real cb = x ? b.bounds.center().x : b.bounds.center().y;
                             return ca != cb ? ca < cb :
                                    (a.proxy != b.proxy ? a.proxy < b.proxy : a.node < b.node);
                         });
    const i32 pair = allocate_pair(new_tree);
    build_subtree(old_tree, items, first, middle, new_tree, pair, target);
    build_subtree(old_tree, items, middle, last, new_tree, pair + 1, target);
    new_tree.nodes[target] = internal_node(new_tree, pair);
    new_tree.parents[target] = parent;
}
void DynamicAabbTree::rebuild_tree(u8 type, bool full) {
    const Tree& old_tree = trees_[type];
    if (!full && (old_tree.nodes.empty() || is_empty(old_tree.nodes[0]) ||
                  (!moved(old_tree.nodes[0]) && old_tree.ordered))) return;
    auto& items = rebuild_items_;
    items.clear();
    if (full) {
        for (i32 i = 0; static_cast<std::size_t>(i) < proxies_.size(); ++i) {
            const auto& p = proxies_[i];
            if (p.live && p.tree == type)
                items.push_back({expanded_bounds(p.exact, p.flags, proxy_extension_ratio_),
                                 -1, i});
        }
    } else {
        gather_build_items(old_tree, 0, items, false);
    }
    Tree& next = rebuild_scratch_;
    next.nodes.assign(2, empty_node());
    next.parents.assign(2, -1);
    next.free_pairs.clear();
    next.ordered = true;
    if (!items.empty())
        build_subtree(old_tree, items, 0, items.size(), next, 0, -1);
    std::swap(trees_[type], next);
}
void DynamicAabbTree::rebuild() {
    for (u8 type = 0; type < 3; ++type) rebuild_tree(type, true);
}
void DynamicAabbTree::rebuild_changed() {
    rebuild_tree(1, false);
    rebuild_tree(2, false);
}
bool DynamicAabbTree::inspect(std::array<QualityStats, 3>& quality) const {
    quality = {};
    for (u8 type = 0; type < 3; ++type) {
        const Tree& tree = trees_[type];
        auto& result = quality[type];
        u32 expected_leaves = 0;
        for (const auto& proxy : proxies_)
            expected_leaves += proxy.live && proxy.node >= 0 && proxy.tree == type ? 1u : 0u;
        if (tree.nodes.empty() || is_empty(tree.nodes[0])) {
            if (expected_leaves != 0) return false;
            continue;
        }
        if (tree.nodes.size() != tree.parents.size() || tree.parents[0] != -1)
            return false;
        std::vector<u8> seen(tree.nodes.size(), 0);
        std::vector<i32> order;
        std::vector<std::pair<i32, i32>> stack{{0, -1}};
        double internal_perimeter = 0, overlap_area = 0, parent_area = 0;
        while (!stack.empty()) {
            const auto [index, parent] = stack.back();
            stack.pop_back();
            if (index < 0 || static_cast<std::size_t>(index) >= tree.nodes.size() ||
                seen[index] || tree.parents[index] != parent)
                return false;
            seen[index] = 1;
            order.push_back(index);
            ++result.reachable_nodes;
            const Node& node = tree.nodes[index];
            if (is_empty(node)) return false;
            if (is_leaf(node)) {
                const i32 proxy = proxy_index(node);
                if (proxy < 0 || static_cast<std::size_t>(proxy) >= proxies_.size() ||
                    !proxies_[proxy].live || proxies_[proxy].tree != type ||
                    proxies_[proxy].node != index ||
                    !contains(node.bounds, proxies_[proxy].exact))
                    return false;
                ++result.leaf_count;
                continue;
            }
            const i32 pair = child(node);
            if (pair < 2 || (pair & 1) != 0 ||
                static_cast<std::size_t>(pair + 1) >= tree.nodes.size())
                return false;
            if (index != 0) internal_perimeter += perimeter(node.bounds);
            const Aabb& a = tree.nodes[pair].bounds;
            const Aabb& b = tree.nodes[pair + 1].bounds;
            if (node.bounds != Aabb::merge(a, b)) return false;
            const double width = std::max(real(0), std::min(a.max.x, b.max.x) -
                                                   std::max(a.min.x, b.min.x));
            const double height = std::max(real(0), std::min(a.max.y, b.max.y) -
                                                    std::max(a.min.y, b.min.y));
            overlap_area += width * height;
            parent_area += static_cast<double>(node.bounds.max.x - node.bounds.min.x) *
                           static_cast<double>(node.bounds.max.y - node.bounds.min.y);
            stack.emplace_back(pair, index);
            stack.emplace_back(pair + 1, index);
        }
        if (result.leaf_count != expected_leaves ||
            result.reachable_nodes != 2 * result.leaf_count - 1)
            return false;
        std::vector<i32> heights(tree.nodes.size(), 0);
        for (auto it = order.rbegin(); it != order.rend(); ++it) {
            const i32 index = *it;
            const Node& node = tree.nodes[index];
            if (!is_leaf(node)) {
                const i32 pair = child(node);
                heights[index] = 1 + std::max(heights[pair], heights[pair + 1]);
                if (node.height_or_shape != heights[index]) return false;
            }
        }
        result.height = static_cast<u32>(heights[0]);
        const double root_perimeter = perimeter(tree.nodes[0].bounds);
        result.area_ratio = root_perimeter > 0 ? internal_perimeter / root_perimeter : 0;
        result.sibling_overlap_ratio = parent_area > 0 ? overlap_area / parent_area : 0;
    }
    return true;
}
void DynamicAabbTree::query_tree(u8 type, u32 id, const Aabb& bounds, u32 flags,
                                  std::vector<BroadPhase::Pair>& out,
                                  std::vector<i32>* scratch,
                                  BroadPhase::QueryStats* stats,
                                  std::span<const u32> marks, u32 epoch) const {
    const Tree& tree = trees_[type];
    if (tree.nodes.empty() || is_empty(tree.nodes[0])) return;
    std::vector<i32> local;
    auto& stack = scratch ? *scratch : local;
    stack.clear();
    stack.push_back(0);
    if (stats) ++stats->nodes_visited;
    while (!stack.empty()) {
        const i32 index = stack.back();
        stack.pop_back();
        const Node& node = tree.nodes[index];
        if (!Aabb::overlap(node.bounds, bounds)) continue;
        if (!is_leaf(node)) {
            const i32 pair = child(node);
            if (stats) stats->nodes_visited += 2;
            if (Aabb::overlap(tree.nodes[pair].bounds, bounds)) stack.push_back(pair);
            if (Aabb::overlap(tree.nodes[pair + 1].bounds, bounds)) stack.push_back(pair + 1);
            continue;
        }
        if (stats) ++stats->leaves_visited;
        const auto& proxy = proxies_[proxy_index(node)];
        if (id == proxy.id || !should_pair(flags, proxy.flags) ||
            (proxy.id < marks.size() && marks[proxy.id] == epoch && proxy.id < id) ||
            !Aabb::overlap(bounds, proxy.exact))
            continue;
        out.push_back(id < proxy.id ? BroadPhase::Pair{id, proxy.id}
                                    : BroadPhase::Pair{proxy.id, id});
        if (stats) ++stats->candidates_emitted;
    }
}
void DynamicAabbTree::query_pairs(u32 id, const Aabb& bounds, u32 flags,
                                   std::vector<BroadPhase::Pair>& out,
                                   std::vector<i32>* scratch,
                                   BroadPhase::QueryStats* stats,
                                   std::span<const u32> marks, u32 epoch) const {
    if (!bounds.is_valid()) return;
    for (u8 type = 0; type < 3; ++type)
        if (type != 0 || !(flags & BroadPhase::PROXY_STATIC))
            query_tree(type, id, bounds, flags, out, scratch, stats, marks, epoch);
}
bool DynamicAabbTree::query_proxy_pairs(i32 handle, std::vector<BroadPhase::Pair>& out,
                                         std::vector<i32>* scratch,
                                         BroadPhase::QueryStats* stats,
                                         std::span<const u32> marks, u32 epoch) const {
    if (!valid_handle(handle) || proxies_[handle].node < 0) return false;
    const auto& p = proxies_[handle];
    query_pairs(p.id, p.exact, p.flags, out, scratch, stats, marks, epoch);
    return true;
}
void DynamicAabbTree::build_pairs(std::vector<BroadPhase::Pair>& out,
                                   std::vector<i32>* scratch, u64* ordering_sorts) const {
    out.clear();
    for (const auto& p : proxies_)
        if (p.live && p.node >= 0)
            query_pairs(p.id, p.exact, p.flags, out, scratch);
    if (!std::is_sorted(out.begin(), out.end())) {
        std::sort(out.begin(), out.end());
        if (ordering_sorts) ++*ordering_sorts;
    }
    out.erase(std::unique(out.begin(), out.end()), out.end());
}
void DynamicAabbTree::query_aabb(const Aabb& bounds, std::vector<u32>& out) const {
    out.clear();
    if (!bounds.is_valid()) return;
    for (const Tree& tree : trees_) {
        if (tree.nodes.empty() || is_empty(tree.nodes[0])) continue;
        std::vector<i32> stack{0};
        while (!stack.empty()) {
            const i32 index = stack.back();
            stack.pop_back();
            const Node& node = tree.nodes[index];
            if (!Aabb::overlap(node.bounds, bounds)) continue;
            if (is_leaf(node)) {
                const auto& p = proxies_[proxy_index(node)];
                if (Aabb::overlap(p.exact, bounds)) out.push_back(p.id);
            } else {
                const i32 pair = child(node);
                stack.push_back(pair);
                stack.push_back(pair + 1);
            }
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
}
BroadPhase::RayCastStats DynamicAabbTree::ray_cast(const BroadPhase::RayCastInput& in,
                                                    BroadPhase::RayCastCallback callback,
                                                    void* context,
                                                    std::vector<i32>* scratch) const {
    BroadPhase::RayCastStats stats;
    const real length2 = math::dot(in.translation, in.translation);
    if (!callback || !std::isfinite(in.origin.x) || !std::isfinite(in.origin.y) ||
        !std::isfinite(in.translation.x) || !std::isfinite(in.translation.y) ||
        !std::isfinite(in.max_fraction) || in.max_fraction <= 0 ||
        length2 <= config::EPSILON * config::EPSILON) return stats;
    if (scratch) scratch->clear();
    const Vec2 direction = in.translation / std::sqrt(length2);
    const Vec2 normal{-direction.y, direction.x};
    const Vec2 absolute{std::abs(normal.x), std::abs(normal.y)};
    real fraction = in.max_fraction;
    const auto segment = [&] {
        const Vec2 end = in.origin + in.translation * fraction;
        return Aabb{{std::min(in.origin.x, end.x), std::min(in.origin.y, end.y)},
                    {std::max(in.origin.x, end.x), std::max(in.origin.y, end.y)}};
    };
    Aabb box = segment();
    std::vector<std::pair<u8, i32>> stack;
    for (u8 type = 0; type < 3; ++type)
        if (!trees_[type].nodes.empty() && !is_empty(trees_[type].nodes[0]))
            stack.emplace_back(type, 0);
    BroadPhase::RayCastInput clipped = in;
    while (!stack.empty()) {
        const auto [type, index] = stack.back();
        stack.pop_back();
        const Tree& tree = trees_[type];
        const Node& node = tree.nodes[index];
        ++stats.nodes_visited;
        if (!Aabb::overlap(node.bounds, box)) continue;
        const Vec2 center = node.bounds.center();
        const Vec2 half = node.bounds.extent();
        const Vec2 delta = in.origin - center;
        if (std::abs(math::dot(normal, delta)) > math::dot(absolute, half)) continue;
        if (is_leaf(node)) {
            const auto& proxy = proxies_[proxy_index(node)];
            if (!Aabb::overlap(proxy.exact, box)) continue;
            const Vec2 exact_delta = in.origin - proxy.exact.center();
            if (std::abs(math::dot(normal, exact_delta)) >
                math::dot(absolute, proxy.exact.extent())) continue;
            clipped.max_fraction = fraction;
            const real result = callback(clipped, proxy.id, context);
            ++stats.leaves_visited;
            if (result == 0) return stats;
            if (std::isfinite(result) && result > 0 && result <= fraction) {
                fraction = result;
                box = segment();
            }
        } else {
            const i32 pair = child(node);
            const real left = distance_squared(tree.nodes[pair].bounds.center(), in.origin);
            const real right = distance_squared(tree.nodes[pair + 1].bounds.center(), in.origin);
            if (left < right) {
                stack.emplace_back(type, pair + 1);
                stack.emplace_back(type, pair);
            } else {
                stack.emplace_back(type, pair);
                stack.emplace_back(type, pair + 1);
            }
        }
    }
    return stats;
}
void DynamicAabbTree::append_state_hash(detail::StateHash& hash) const {
    hash.add(proxy_extension_ratio_);
    hash.add(static_cast<u64>(proxy_count_));
    hash.add(static_cast<u64>(proxies_.size()));
    for (const auto& proxy : proxies_) {
        hash.add(proxy.live);
        if (proxy.live) {
            hash.add(proxy.exact); hash.add(proxy.id); hash.add(proxy.flags);
            hash.add(proxy.node); hash.add(proxy.tree);
        }
    }
    for (const auto& tree : trees_) {
        hash.add(tree.ordered);
        hash.add(static_cast<u64>(tree.nodes.size()));
        for (std::size_t i = 0; i < tree.nodes.size(); ++i) {
            const auto& node = tree.nodes[i];
            hash.add(is_empty(node));
            if (!is_empty(node)) {
                hash.add(node.bounds); hash.add(node.flag_index);
                hash.add(node.height_or_shape); hash.add(tree.parents[i]);
            }
        }
        hash.sequence(tree.free_pairs);
    }
}

} // namespace sqp
