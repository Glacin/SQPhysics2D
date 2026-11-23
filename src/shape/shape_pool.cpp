// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <array>
#include <limits>
#include <utility>
#include <vector>

#include "sqphysics2d/detail/data_storage.h"
#include "world/state_hash.h"

namespace sqp::detail {

enum class ShapeColumn : u8 {
    CenterX,
    CenterY,
    Area,
    InertiaPerUnitMass,
    BoundRadius,
    InnerBoundRadius,
    Density,
    CircleRadius,
    RestitutionImpactScale,
    Count
};

inline constexpr u32 invalid_shape_dense_index = std::numeric_limits<u32>::max();
inline constexpr std::size_t shape_field_count = static_cast<std::size_t>(ShapeColumn::Count);
using ShapeValues = std::array<real, shape_field_count>;

struct ShapeSparseSlot {
    u32 dense_index = invalid_shape_dense_index;
    u32 generation = 1;

    void release() {
        dense_index = invalid_shape_dense_index;
        if (++generation == 0)
            generation = 1;
    }
};

// Both shape kinds use the same columns; polygon rows leave CircleRadius unused.
struct ShapeStorage {
    std::array<std::vector<real>, shape_field_count> columns;
    std::vector<u32> material_indices;
    std::vector<u32> dense_to_sparse;
    std::vector<ShapeSparseSlot> sparse;
    std::vector<u32> free_list;

    u32 find(ShapeHandle handle) const {
        if (handle.index >= sparse.size())
            return invalid_shape_dense_index;
        const auto& slot = sparse[handle.index];
        return slot.generation == handle.generation ? slot.dense_index : invalid_shape_dense_index;
    }
    u32 size() const {
        return static_cast<u32>(dense_to_sparse.size());
    }
    void reserve(u32 capacity);
    ShapeHandle insert(ShapeKind kind, const ShapeValues& values, u32 material_index);
    void write(u32 dense, const ShapeValues& values, u32 material_index);
    void erase(ShapeHandle handle, u32 dense);
    void clear();
};

struct PolygonStorage : ShapeStorage {
    std::vector<ConvexPoly> geometry;
};

void ShapeStorage::reserve(u32 capacity) {
    for (auto& column : columns)
        column.reserve(capacity);
    material_indices.reserve(capacity);
    dense_to_sparse.reserve(capacity);
    sparse.reserve(capacity);
    free_list.reserve(capacity);
}

ShapeHandle ShapeStorage::insert(ShapeKind kind, const ShapeValues& values, u32 material_index) {
    if (dense_to_sparse.size() >= invalid_shape_dense_index)
        return {};
    u32 index;
    if (free_list.empty()) {
        if (sparse.size() >= invalid_shape_dense_index)
            return {};
        index = static_cast<u32>(sparse.size());
        sparse.push_back({});
    } else {
        index = free_list.back();
        free_list.pop_back();
    }
    const u32 dense = size();
    for (std::size_t field = 0; field < columns.size(); ++field)
        columns[field].push_back(values[field]);
    material_indices.push_back(material_index);
    dense_to_sparse.push_back(index);
    sparse[index].dense_index = dense;
    return {index, sparse[index].generation, kind};
}

void ShapeStorage::write(u32 dense, const ShapeValues& values, u32 material_index) {
    for (std::size_t field = 0; field < columns.size(); ++field)
        columns[field][dense] = values[field];
    material_indices[dense] = material_index;
}

void ShapeStorage::erase(ShapeHandle handle, u32 dense) {
    const u32 last = size() - 1;
    if (dense != last) {
        for (auto& column : columns)
            column[dense] = column[last];
        material_indices[dense] = material_indices[last];
        const u32 moved = dense_to_sparse[last];
        dense_to_sparse[dense] = moved;
        sparse[moved].dense_index = dense;
    }
    for (auto& column : columns)
        column.pop_back();
    material_indices.pop_back();
    dense_to_sparse.pop_back();
    sparse[handle.index].release();
    free_list.push_back(handle.index);
}

void ShapeStorage::clear() {
    for (auto& column : columns)
        column.clear();
    material_indices.clear();
    for (const u32 index : dense_to_sparse)
        sparse[index].release();
    dense_to_sparse.clear();
    free_list.clear();
    free_list.reserve(sparse.size());
    for (u32 i = 0; i < static_cast<u32>(sparse.size()); ++i)
        free_list.push_back(i);
}

} // namespace sqp::detail

namespace sqp {

struct ShapePool::Impl {
    MaterialPool materials;
    detail::ShapeStorage circles;
    detail::PolygonStorage polygons;

    const detail::ShapeStorage* resolve(ShapeHandle handle, u32& dense) const {
        const detail::ShapeStorage* storage;
        switch (handle.kind) {
        case ShapeKind::Circle:
            storage = &circles;
            break;
        case ShapeKind::Polygon:
            storage = &polygons;
            break;
        default:
            return nullptr;
        }
        dense = storage->find(handle);
        return dense == detail::invalid_shape_dense_index ? nullptr : storage;
    }
    detail::ShapeStorage* resolve(ShapeHandle handle, u32& dense) {
        return const_cast<detail::ShapeStorage*>(
            static_cast<const Impl&>(*this).resolve(handle, dense));
    }
};

namespace {
using detail::ShapeColumn;

PolygonVertices transformed_vertices(const PolygonShapeDef& def) {
    PolygonVertices result = def.vertices;
    if (def.vertices.overflowed())
        return result;
    const real cosine = math::cos(def.local_rotation);
    const real sine = math::sin(def.local_rotation);
    for (Vec2& vertex : result)
        vertex = math::rotate_vector(vertex, cosine, sine) + def.local_offset;
    return result;
}

struct GeometryProperties {
    real area = 0, inertia = 0, bound = 0, inner_bound = 0;
};

GeometryProperties geometry_properties(const Circle& circle) {
    if (circle.radius <= real(0))
        return {};
    return {math::PI * circle.radius * circle.radius, real(0.5) * circle.radius * circle.radius,
            circle.radius, circle.radius};
}

GeometryProperties geometry_properties(const ConvexPoly& polygon) {
    GeometryProperties result;
    if (polygon.local.size() < 3)
        return result;
    real twice_area = 0, inertia_numerator = 0, radius_squared = 0;
    result.inner_bound = std::numeric_limits<real>::max();
    for (std::size_t i = 0; i < polygon.local.size(); ++i) {
        const Vec2& a = polygon.local[i];
        const Vec2& b = polygon.local[(i + 1) % polygon.local.size()];
        const Vec2 edge = b - a;
        const real edge_length = edge.length();
        if (edge_length <= config::EPSILON)
            return {};
        const real cross = math::cross(a, b);
        twice_area += cross;
        inertia_numerator += cross * (math::dot(a, a) + math::dot(a, b) + math::dot(b, b));
        radius_squared = std::max(radius_squared, a.length_squared());
        result.inner_bound = std::min(result.inner_bound, math::cross(edge, -a) / edge_length);
    }
    if (twice_area <= config::EPSILON)
        return {};
    result.area = twice_area * real(0.5);
    result.inertia = inertia_numerator / (real(6) * twice_area);
    result.bound = math::sqrt(radius_squared);
    result.inner_bound = std::max(real(0), result.inner_bound);
    return result;
}

// Use the geometry formulas once per successful creation/replacement.
template <class Geometry>
detail::ShapeValues cached_values(const Geometry& geometry, real density,
                                  real restitution_impact_scale, real radius = real(0)) {
    const GeometryProperties computed = geometry_properties(geometry);
    detail::ShapeValues values{};
    const auto put = [&](ShapeColumn field, real value) {
        values[static_cast<std::size_t>(field)] = value;
    };
    const Vec2 center = [&]() -> Vec2 {
        if constexpr (std::is_same_v<Geometry, Circle>)
            return geometry.center;
        else
            return geometry.center_offset();
    }();
    put(ShapeColumn::CenterX, center.x);
    put(ShapeColumn::CenterY, center.y);
    put(ShapeColumn::Area, computed.area);
    put(ShapeColumn::InertiaPerUnitMass, computed.inertia);
    put(ShapeColumn::BoundRadius, computed.bound);
    put(ShapeColumn::InnerBoundRadius, computed.inner_bound);
    put(ShapeColumn::Density, std::max(real(0), density));
    put(ShapeColumn::CircleRadius, radius);
    put(ShapeColumn::RestitutionImpactScale,
        std::max(real(1) - config::EPSILON, restitution_impact_scale));
    return values;
}
} // namespace

ShapePool::ShapePool() : impl_(std::make_unique<Impl>()) {}
ShapePool::~ShapePool() = default;
ShapePool::ShapePool(const ShapePool& other) : impl_(std::make_unique<Impl>(*other.impl_)) {}
ShapePool& ShapePool::operator=(const ShapePool& other) {
    if (this != &other)
        impl_ = std::make_unique<Impl>(*other.impl_);
    return *this;
}
ShapePool::ShapePool(ShapePool&&) noexcept = default;
ShapePool& ShapePool::operator=(ShapePool&&) noexcept = default;

void ShapePool::clear() {
    impl_->circles.clear();
    impl_->polygons.clear();
    impl_->polygons.geometry.clear();
    impl_->materials.reset();
}

void ShapePool::reserve(u32 circle_capacity, u32 polygon_capacity) {
    impl_->circles.reserve(circle_capacity);
    impl_->polygons.reserve(polygon_capacity);
    impl_->polygons.geometry.reserve(polygon_capacity);
}

ShapeHandle ShapePool::add_circle(const CircleShapeDef& def) {
    const Circle circle(def.local_offset, def.radius);
    if (circle.radius <= real(0) || !impl_->materials.contains(def.material_index))
        return {};
    const auto values =
        cached_values(circle, def.density, def.restitution_impact_scale, circle.radius);
    return impl_->circles.insert(ShapeKind::Circle, values, def.material_index);
}

ShapeHandle ShapePool::add_polygon(const PolygonShapeDef& def) {
    const ConvexPoly polygon(transformed_vertices(def));
    if (polygon.local.size() < 3 || !impl_->materials.contains(def.material_index))
        return {};
    const auto values = cached_values(polygon, def.density, def.restitution_impact_scale);
    const ShapeHandle handle =
        impl_->polygons.insert(ShapeKind::Polygon, values, def.material_index);
    if (sqp::is_valid(handle))
        impl_->polygons.geometry.push_back(polygon);
    return handle;
}

bool ShapePool::get_circle(ShapeHandle handle, Circle& out) const {
    u32 dense;
    const auto* storage = static_cast<const Impl&>(*impl_).resolve(handle, dense);
    if (!storage || handle.kind != ShapeKind::Circle)
        return false;
    Circle value;
    value.center = {storage->columns[static_cast<std::size_t>(ShapeColumn::CenterX)][dense],
                    storage->columns[static_cast<std::size_t>(ShapeColumn::CenterY)][dense]};
    value.radius = storage->columns[static_cast<std::size_t>(ShapeColumn::CircleRadius)][dense];
    out = value;
    return true;
}

bool ShapePool::get_polygon(ShapeHandle handle, PolygonVertices& out) const {
    const auto* polygon = get_convex_polygon(handle);
    if (!polygon)
        return false;
    PolygonVertices result;
    for (const Vec2& vertex : polygon->local)
        result.push_back(polygon->center_offset() + vertex);
    out = result;
    return true;
}

bool ShapePool::get_properties(ShapeHandle handle, ShapeProperties& out) const {
    u32 dense;
    const auto* storage = static_cast<const Impl&>(*impl_).resolve(handle, dense);
    if (!storage)
        return false;
    const auto value = [&](ShapeColumn field) {
        return storage->columns[static_cast<std::size_t>(field)][dense];
    };
    out = {{value(ShapeColumn::CenterX), value(ShapeColumn::CenterY)},
           value(ShapeColumn::Area),
           value(ShapeColumn::InertiaPerUnitMass),
           value(ShapeColumn::BoundRadius),
           value(ShapeColumn::InnerBoundRadius),
           value(ShapeColumn::Density),
           value(ShapeColumn::RestitutionImpactScale),
           storage->material_indices[dense]};
    return true;
}

bool ShapePool::set_density(ShapeHandle handle, real density) {
    u32 dense;
    auto* storage = impl_->resolve(handle, dense);
    if (!storage)
        return false;
    storage->columns[static_cast<std::size_t>(ShapeColumn::Density)][dense] =
        std::max(real(0), density);
    return true;
}

bool ShapePool::set_restitution_impact_scale(ShapeHandle handle, real restitution_impact_scale) {
    u32 dense;
    auto* storage = impl_->resolve(handle, dense);
    if (!storage)
        return false;
    storage->columns[static_cast<std::size_t>(ShapeColumn::RestitutionImpactScale)][dense] =
        std::max(real(1) - config::EPSILON, restitution_impact_scale);
    return true;
}

bool ShapePool::get_material_index(ShapeHandle handle, u32& out) const {
    u32 dense;
    const auto* storage = static_cast<const Impl&>(*impl_).resolve(handle, dense);
    if (!storage)
        return false;
    out = storage->material_indices[dense];
    return true;
}

bool ShapePool::get_material(ShapeHandle handle, PhysicsMaterial& out) const {
    u32 dense;
    const auto* storage = static_cast<const Impl&>(*impl_).resolve(handle, dense);
    if (!storage)
        return false;
    return impl_->materials.get(storage->material_indices[dense], out);
}

bool ShapePool::set_material_index(ShapeHandle handle, u32 index) {
    u32 dense;
    auto* storage = impl_->resolve(handle, dense);
    if (!storage || !impl_->materials.contains(index))
        return false;
    storage->material_indices[dense] = index;
    return true;
}

MaterialPool& ShapePool::materials() {
    return impl_->materials;
}

const MaterialPool& ShapePool::materials() const {
    return impl_->materials;
}

bool ShapePool::update_circle(ShapeHandle handle, const CircleShapeDef& def) {
    if (handle.kind != ShapeKind::Circle)
        return false;
    u32 dense;
    auto* storage = impl_->resolve(handle, dense);
    const Circle circle(def.local_offset, def.radius);
    if (!storage || circle.radius <= real(0) || !impl_->materials.contains(def.material_index))
        return false;
    const auto values =
        cached_values(circle, def.density, def.restitution_impact_scale, circle.radius);
    storage->write(dense, values, def.material_index);
    return true;
}

bool ShapePool::update_polygon(ShapeHandle handle, const PolygonShapeDef& def) {
    if (handle.kind != ShapeKind::Polygon)
        return false;
    u32 dense;
    auto* storage = impl_->resolve(handle, dense);
    if (!storage)
        return false;
    const ConvexPoly polygon(transformed_vertices(def));
    if (polygon.local.size() < 3 || !impl_->materials.contains(def.material_index))
        return false;
    const auto values = cached_values(polygon, def.density, def.restitution_impact_scale);
    impl_->polygons.geometry[dense] = polygon;
    storage->write(dense, values, def.material_index);
    return true;
}

bool ShapePool::remove(ShapeHandle handle) {
    u32 dense;
    auto* storage = impl_->resolve(handle, dense);
    if (!storage)
        return false;
    if (handle.kind == ShapeKind::Polygon) {
        auto& geometry = impl_->polygons.geometry;
        if (dense + 1 != geometry.size())
            geometry[dense] = std::move(geometry.back());
        geometry.pop_back();
    }
    storage->erase(handle, dense);
    return true;
}

bool ShapePool::is_valid(ShapeHandle handle) const {
    u32 dense;
    return static_cast<const Impl&>(*impl_).resolve(handle, dense) != nullptr;
}

ShapeKind ShapePool::kind(ShapeHandle handle) const {
    return is_valid(handle) ? handle.kind : ShapeKind::None;
}

u32 ShapePool::shape_count() const {
    const std::size_t count =
        static_cast<std::size_t>(impl_->circles.size()) + impl_->polygons.size();
    return count >= INVALID_SHAPE_INDEX ? INVALID_SHAPE_INDEX : static_cast<u32>(count);
}

const ConvexPoly* ShapePool::get_convex_polygon(ShapeHandle handle) const {
    if (handle.kind != ShapeKind::Polygon)
        return nullptr;
    u32 dense;
    if (!static_cast<const Impl&>(*impl_).resolve(handle, dense))
        return nullptr;
    return &impl_->polygons.geometry[dense];
}

void ShapePool::append_state_hash(detail::StateHash& hash) const {
    const auto& materials = impl_->materials;
    hash.add(materials.size());
    for (u32 i = 0; i < materials.size(); ++i) {
        PhysicsMaterial material;
        materials.get(i, material);
        hash.add(material);
    }
    const auto storage = [&](const detail::ShapeStorage& value) {
        for (const auto& column : value.columns)
            hash.sequence(column);
        hash.sequence(value.material_indices);
        hash.sequence(value.dense_to_sparse);
        hash.add(static_cast<u64>(value.sparse.size()));
        for (const auto& slot : value.sparse) {
            hash.add(slot.dense_index); hash.add(slot.generation);
        }
        hash.sequence(value.free_list);
    };
    storage(impl_->circles);
    storage(impl_->polygons);
    hash.add(static_cast<u64>(impl_->polygons.geometry.size()));
    for (const auto& polygon : impl_->polygons.geometry) {
        hash.add(polygon.center_offset());
        hash.sequence(polygon.local);
    }
}

} // namespace sqp
