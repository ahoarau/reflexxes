// Copyright 2026 Intrinsic Innovation LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef INTRINSIC_ICON_REFLEXXES_INTERNAL_MATH_H_
#define INTRINSIC_ICON_REFLEXXES_INTERNAL_MATH_H_

#include <math.h>

#include <algorithm>

#include "intrinsic/icon/reflexxes/internal/constants.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

// Positive zero value to use to prevent numerical errors.
constexpr double kPositiveZero = 1.0e-50;

// Sign function (integer)
inline int GetSign(const double A) {
  if (A < 0.0) {
    return -1;
  } else {
    return 1;
  }
}

// Calculates the real square root of a given value. If the value is
// negative a value of almost zero will be returned.
inline double GetSqrt(const double value) {
  return ((value <= 0.0) ? (kPositiveZero) : (sqrt(value)));
}

// Returns true if the values are within epsilon of each other.
inline bool EpsilonEqual(const double a, const double b, const double epsilon) {
  return fabs(a - b) <= epsilon;
}

// Returns true if the value is small enough to be considered zero.
inline double GetNearZero(const double value) {
  return static_cast<double>(fabs(value) <= kPositiveZero);
}

// Returns x^2
inline double Power2(const double x) { return (x * x); }

// Returns x^3
inline double Power3(const double x) { return (x * x * x); }

// Returns x^4
inline double Power4(const double x) { return (x * x * x * x); }

// The min/max range type
struct Range {
  double min;
  double max;
};

// Clips x to [limits.first, limits.second]
inline double ClipToRange(const double x, const Range& limits) {
  if (EpsilonEqual(limits.min, limits.max, kGeneralEqualityEpsilon))
    return limits.max;
  return std::clamp(x, limits.min, limits.max);
}

// Clips x within [min, max]
inline double ClipToRange(const double x, const double min, const double max) {
  if (EpsilonEqual(min, max, kGeneralEqualityEpsilon)) return max;
  return std::clamp(x, min, max);
}

inline Range ClipToRange(const Range& x, const Range& limits) {
  return {ClipToRange(x.min, limits), ClipToRange(x.max, limits)};
}

inline Range ClipToRange(const Range& x, const double min, const double max) {
  return ClipToRange(x, {min, max});
}

// Returns the extent of the range.
inline double GetRangeSize(const Range& x) { return x.max - x.min; }

// Returns true if the given value is between [min, max]
inline bool IsInRange(const Range& range, const double x) {
  return x >= range.min && x <= range.max;
}

// If the given range has a range larger than a relative epsilon, shrink the
// range by that epsilon.  This is useful when the range is a limit on some
// boundary condition, and a small change in value will overcome that boundary.
inline Range ShrinkRange(const Range& value_min_max) {
  // shrink the range of values
  auto [value_min, value_max] = value_min_max;
  double offset_value =
      kRelativeLimitOffset * (value_max - value_min) + kAbsoluteLimitOffset;
  if ((value_max - value_min) >= 2 * offset_value) {
    return {value_min + offset_value, value_max - offset_value};
  }
  return value_min_max;
}

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_INTERNAL_MATH_H_
