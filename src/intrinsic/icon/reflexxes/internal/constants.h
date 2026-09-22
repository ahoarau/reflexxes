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

#ifndef INTRINSIC_ICON_REFLEXXES_INTERNAL_CONSTANTS_H_
#define INTRINSIC_ICON_REFLEXXES_INTERNAL_CONSTANTS_H_

#include "intrinsic/icon/reflexxes/constants_external.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

// A value for infinity.
inline constexpr double kInfinity = 1.0e100;

// Stop condition for the Anderson-Bjoerck-King method.
inline constexpr double kABKEpsilon = 1.0e-12;

// Value for the Step2 position error functions to prevent numerical faults.
inline constexpr double kABKFunctionEpsilon = 1.0e-14;

// Epsilon for time offsets in Step 1A and Step 1C.
inline constexpr double kStep1TimeEpsilon = 1.0e-10;

// Relative offset to expand the limit parameters for the Anderson-Bjoerck-King
// method.
inline constexpr double kRelativeLimitOffset = 1.0e-9;

// Absolute offset to expand the limit parameters for the Anderson-Bjoerck-King
// method.
inline constexpr double kAbsoluteLimitOffset = 1.0e-8;

// Epsilon value to consider a value as a zero in order to determine whether
// a motion profile can provide a valid trajectory.
inline constexpr double kValidSolutionEpsilon = 1.0e-3;

// Absolute epsilon to check whether all required vectors
// are collinear
inline constexpr double kAbsPhaseSyncEpsilon = 1.0e-6;

// Relative epsilon to check whether all required vectors
// are collinear.
inline constexpr double kRelPhaseSyncEpsilon = 1.0e-3;

// Absolute epsilon value to increase the position check in case
// of target velocities near the maximum velocity.
inline constexpr double kStep2VMaxEpsilon = 1.0e-9;

// Absolute epsilon to check whether the target position
// has been already reached.
inline constexpr double kAbsoluteStep1B1Epsilon = 1.0e-9;

// Relative epsilon to check whether the target position
// has been already reached.
inline constexpr double kRelativeStep1B1Epsilon = 1.0e-9;

// A relative value, from when on a position error during a root
// calculation in Step 1 shall be considered as an error.
inline constexpr double kRelStep1PositionErrorTolerance = 1.0e-2;

// An absolute value, from when on a position error during a root
// calculation in Step 1 shall be considered as an error.
inline constexpr double kAbsStep1PositionErrorTolerance = 1.0e-1;

// A relative value, from when on a position error during a root
// calculation in Step 2 shall be considered as an error.
inline constexpr double kRelPositionErrorTolerance = 1.0e-2;

// An absolute value, from when on a position error during a root
// calculation in Step 2 shall be considered as an error.
inline constexpr double kAbsPositionErrorTolerance = 1.0e-1;

// A relative value, from when on a velocity error during a root
// calculation in Step 2 shall be considered as an error.
inline constexpr double kRelVelocityErrorTolerance = 1.0e-2;

// An absolute value, from when on a velocity error during a root
// calculation in Step 2 shall be considered as an error.
inline constexpr double kAbsVelocityErrorTolerance = 1.0e-1;

// An absolute value, from when on an acceleration error during a root
// calculation in Step 2 shall be considered as an error.
inline constexpr double kAbsAccelerationErrorTolerance = 1.0e-1;

// An absolute value, from when on a time error during a root
// calculation in Step 2 shall be considered as an error.
inline constexpr double kAbsTimeErrorTolerance = 1.0e-1;

// Absolute epsilon value for time to be considered equal to zero.
inline constexpr double kAbsZeroTimeEpsilon = 1.0e-8;

// The value of minimum synchronization time exceeding which error
// is reported.
inline constexpr double kMaxAllowedSyncTimeSeconds = 1e6;

// Constant values used by the expressions for polynomial coefficients.
inline constexpr double kOneThird = 1.0 / 3.0;
inline constexpr double kOneSixth = 1.0 / 6.0;
inline constexpr double kOneTwentyFourth = 1.0 / 24.0;

// Constant used to expand the parameter limits in order to
// check for valid solutions, if needed.
inline constexpr double kParamLimitOffset = 1.0e-12;

// Constant for state solvers to generally consider something "equal".
inline constexpr double kGeneralEqualityEpsilon = 1.0e-12;

// Threshold for maximum position limit violation to consider the limit
// breached.
inline constexpr double kPositionLimitBreachedThreshold = 1e-12;

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_INTERNAL_CONSTANTS_H_
