// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>

#include "sqphysics2d/detail/types.h"

namespace sqp::config {

// Small-value tolerance used by geometric and numerical helpers.
#if defined(SQP_REAL_DOUBLE)
inline constexpr real EPSILON = real(1e-12);
#else
inline constexpr real EPSILON = real(1e-6);
#endif

// Default collision margin in world units.
inline constexpr real MARGIN = real(0.1);

// Simulation clock precision and nominal tick rate, expressed as powers of two.
inline constexpr u64 TU_PER_SECOND_SHIFT = 24;
inline constexpr u64 TICKS_PER_SECOND_SHIFT = 6;
static_assert(TICKS_PER_SECOND_SHIFT <= TU_PER_SECOND_SHIFT);

// Number of simulation time units in one second and in one nominal tick.
inline constexpr tick_t TU_PER_SECOND = tick_t(1) << TU_PER_SECOND_SHIFT;
inline constexpr tick_t TICKS_PER_SECOND = tick_t(1) << TICKS_PER_SECOND_SHIFT;
inline constexpr tick_t TU_PER_TICK = tick_t(1) << (TU_PER_SECOND_SHIFT - TICKS_PER_SECOND_SHIFT);
static_assert(TU_PER_TICK > 0);

// Number of power-of-two period buckets in the island scheduler.
inline constexpr std::size_t PERIOD_BUCKET_COUNT = 32;
static_assert(PERIOD_BUCKET_COUNT > 0);
static_assert(PERIOD_BUCKET_COUNT > TU_PER_SECOND_SHIFT - TICKS_PER_SECOND_SHIFT);
static_assert(PERIOD_BUCKET_COUNT <= 256); // IslandPool stores the bucket in a u8.

// Maximum number of vertices stored by PolygonVertices.
inline constexpr std::size_t MAX_POLYGON_VERTICES = 8;

// Scalar implementation currently imposes no extra SIMD alignment.
inline constexpr u32 SIMD_ALIGNMENT = 1;

// namespace sqp::config
} // namespace sqp::config
