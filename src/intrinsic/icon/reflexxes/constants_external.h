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

#ifndef INTRINSIC_ICON_REFLEXXES_CONSTANTS_EXTERNAL_H_
#define INTRINSIC_ICON_REFLEXXES_CONSTANTS_EXTERNAL_H_

namespace intrinsic {
namespace reflexxes {

// kMaxDofs defines the maximum number of DOFs the system will accept.
// This constant controls many stack allocated objects, so should be set as low
// as possible within reason.
// For ease of use, the preprocessor define RML_MAX_DOFS can be set from a build
// file to change this value without needing to change the source.
#ifdef RML_MAX_DOFS
inline constexpr int kMaxDofs = RML_MAX_DOFS;
#else
inline constexpr int kMaxDofs = 16;
#endif

// Maximum value for the minimum trajectory execution time.
inline constexpr double kMaxMinExecutionTime = 1e10;

// Lower threshold for kinematic motion constraints used for the scaling.
// the input values.
inline constexpr double kLowerScalingThreshold = 1.0;

// Upper threshold for kinematic motion constraints used for the scaling.
// the input values.
inline constexpr double kUpperScalingThreshold = 100000.0;

// The validity check magnitude range for RMLPosition.
inline constexpr double kPositionValidityMagnitude = 1e8;

// The validity check magnitude range for RMLVelocity.
inline constexpr double kVelocityValidityMagnitude = 1e10;

// The maximum number of polynomials.
inline constexpr unsigned int kMaxNumPolynomials = 14;

// Negative value near zero used for the numerically robust calculation
// square roots.
inline constexpr double kPolynomialRootNegativeZero = 1e-8;

// Positive threshold value to compare current and former input values.
inline constexpr double kInputValueEpsilon = 1e-8;

// Absolute epsilon value to prevent from numerical inaccuracies
// during the check whether the desired target position is exceeded.
inline constexpr double kExceedingTargetPositionEpsilon = 1e-6;

// Absolute maximum positional limit value exceeding which the
// Inputs::IsValid() will return false.
inline constexpr double kMaxPositionalLimit = 1.0e100;

// Absolute maximum jerk limit value exceeding which the Inputs::IsValid() will
// return false.
inline constexpr double kMaxJerkLimit = 1.0e100;

// Determines if profile and decision traces are kept for later debugging.
// Should be off during normal production operation.
inline constexpr bool kEnableTracing = false;

// Size of character arrays used to store decision and profile traces.
inline constexpr unsigned int kTraceSize = (kEnableTracing ? 32 : 0);

}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_CONSTANTS_EXTERNAL_H_
