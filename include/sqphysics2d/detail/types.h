// SPDX-FileCopyrightText: 2026 Glacin
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <type_traits>

namespace sqp {

// Configurable scalar used for simulation calculations.
#if defined(SQP_REAL_DOUBLE)
using real = double;
#else
using real = float;
#endif

// Fixed precision aliases independent of the configured simulation scalar.
using real32 = float;
using real64 = double;
// Signed and unsigned integer aliases with fixed widths.
using i8 = std::int8_t;
using i16 = std::int16_t;
using i32 = std::int32_t;
using i64 = std::int64_t;
using u8 = std::uint8_t;
using u16 = std::uint16_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
// Unsigned simulation-time units; the units per second are defined in build_config.h.
using tick_t = u64;
// Alias for an absolute simulation-time index, still measured in TU.
using tick_index_t = tick_t;

static_assert(std::is_floating_point_v<real>);
static_assert(sizeof(u64) == 8);

// namespace sqp
} // namespace sqp
