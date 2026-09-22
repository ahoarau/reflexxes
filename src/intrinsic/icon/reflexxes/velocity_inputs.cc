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

#include "intrinsic/icon/reflexxes/velocity_inputs.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>

#include "absl/log/check.h"
#include "intrinsic/icon/reflexxes/constants_external.h"
#include "intrinsic/icon/reflexxes/inputs.h"

namespace intrinsic {
namespace reflexxes {

std::pair<Inputs::ErrorCodeForInvalidInputValues, int>
VelocityInputs::CheckForValidity() const {
  if (GetMinimumSynchronizationTime() > kMaxMinExecutionTime) {
    return {ErrorCodeForInvalidInputValues::kMinimumSyncTime, -1};
  }

  for (const DOF& dof : GetDOFs()) {
    if (!dof.selected) {
      continue;
    }

    if ((dof.max_position > kMaxPositionalLimit) ||
        (dof.min_position < -kMaxPositionalLimit)) {
      return {ErrorCodeForInvalidInputValues::kPositionLimitTooLarge,
              dof.index};
    }

    if (dof.max_acceleration <= 0.0) {
      return {ErrorCodeForInvalidInputValues::kMaxAcceleration, dof.index};
    }

    if (dof.min_acceleration >= 0.0) {
      return {ErrorCodeForInvalidInputValues::kMinAcceleration, dof.index};
    }

    if (dof.max_jerk <= 0.0) {
      return {ErrorCodeForInvalidInputValues::kMaxJerk, dof.index};
    }

    if (dof.min_jerk >= 0.0) {
      return {ErrorCodeForInvalidInputValues::kMinJerk, dof.index};
    }

    if ((dof.max_jerk > kMaxJerkLimit) || (dof.min_jerk < -kMaxJerkLimit)) {
      return {ErrorCodeForInvalidInputValues::kJerkLimitTooLarge, dof.index};
    }

    // check for values that can result in numerical instability
    double max_order_acceleration =
        fmax(dof.max_acceleration, fabs(dof.min_acceleration));
    double min_order_acceleration =
        fmin(dof.max_acceleration, fabs(dof.min_acceleration));

    double max_order_jerk = fmax(dof.max_jerk, fabs(dof.min_jerk));
    double min_order_jerk = fmin(dof.max_jerk, fabs(dof.min_jerk));

    // find max order of magnitude
    std::array<double, 6> test{
        max_order_jerk,     max_order_acceleration, fabs(dof.target_velocity),
        fabs(dof.position), fabs(dof.velocity),     fabs(dof.acceleration)};
    double maximum_order_of_magnitude =
        *std::max_element(test.begin(), test.end());

    // And the min order.
    double minimum_order_of_magnitude =
        std::fmin(min_order_acceleration, min_order_jerk);

    // Note: The target velocity value does not have to be checked as we already
    // know that it is less than the maximum velocity value. The alternative
    // target velocity vector does not have to be checked.

    // The value of minimum_order_of_magnitude is greater than
    // zero:
    if ((maximum_order_of_magnitude / minimum_order_of_magnitude) >
        kVelocityValidityMagnitude) {
      return {ErrorCodeForInvalidInputValues::kOrderOfMagnitude, dof.index};
    }
  }

  return {ErrorCodeForInvalidInputValues::kNoError, -1};
}

bool VelocityInputs::CheckValidityOfConstraintsForSingleDOF(
    const Inputs::DOF& dof) {
  return (dof.max_acceleration > 0.0) && (dof.max_jerk > 0.0) &&
         (dof.min_acceleration < 0.0) && (dof.min_jerk < 0.0);
}

bool VelocityInputs::CheckValidityOfConstraints() const {
  for (const DOF& dof : GetDOFs()) {
    if (dof.selected && !CheckValidityOfConstraintsForSingleDOF(dof)) {
      return false;
    }
  }
  return true;
}

void VelocityInputs::DeselectInvalidDofs() {
  for (DOF& dof : GetDOFs()) {
    dof.selected = dof.selected && CheckValidityOfConstraintsForSingleDOF(dof);
  }
}

MaxDOFFixedVector<double> VelocityInputs::ScaleIfNecessary() {
  MaxDOFFixedVector<double> scaling_factors;
  for (DOF& dof : GetDOFs()) {
    double max_order_acceleration =
        fmax(dof.max_acceleration, fabs(dof.min_acceleration));
    double min_order_acceleration =
        fmin(dof.max_acceleration, fabs(dof.min_acceleration));

    double max_order_jerk = fmax(dof.max_jerk, fabs(dof.min_jerk));
    double min_order_jerk = fmin(dof.max_jerk, fabs(dof.min_jerk));

    if (!dof.selected || !CheckValidityOfConstraintsForSingleDOF(dof)) {
      scaling_factors.push_back(1.0);
      continue;
    }

    // First check if jerk or acceleration have exceeded the min scaling
    // threshold.

    if (min_order_jerk < kLowerScalingThreshold) {
      scaling_factors.push_back(kLowerScalingThreshold / min_order_jerk);
      continue;
    }

    if (min_order_acceleration < kLowerScalingThreshold) {
      scaling_factors.push_back(kLowerScalingThreshold /
                                min_order_acceleration);
      continue;
    }

    // Now see if either acceleration or jerk has exceeded the upper scaling
    // threshold.
    double test_scale = 0.0;
    if (max_order_acceleration > kUpperScalingThreshold) {
      test_scale = kUpperScalingThreshold / max_order_acceleration;
    } else {
      if (max_order_jerk > kUpperScalingThreshold) {
        test_scale = kUpperScalingThreshold / max_order_jerk;
      } else {
        scaling_factors.push_back(1.0);
        continue;
      }
    }

    // If the upper thresholds were exceeded, check the lower thresholds again
    // with the new scaling value. If they now exceed, use those instead.
    if ((min_order_jerk * test_scale) < kLowerScalingThreshold) {
      scaling_factors.push_back(kLowerScalingThreshold / min_order_jerk);
      continue;
    }

    if ((min_order_acceleration * test_scale) < kLowerScalingThreshold) {
      scaling_factors.push_back(kLowerScalingThreshold /
                                min_order_acceleration);
      continue;
    }

    scaling_factors.push_back(test_scale);
  }

  // Now scale them all.
  CHECK_EQ(scaling_factors.size(), GetNumberOfDOFs());
  for (Inputs::DOF& dof : GetDOFs()) {
    dof.ScaleValues(scaling_factors[dof.index]);
  }

  return scaling_factors;
}

}  // namespace reflexxes
}  // namespace intrinsic
