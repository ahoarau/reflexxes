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

#include "intrinsic/icon/reflexxes/internal/input_output_utility_functions.h"

#include "absl/log/check.h"
#include "intrinsic/icon/reflexxes/constants_external.h"
#include "intrinsic/icon/reflexxes/inputs.h"
#include "intrinsic/icon/reflexxes/internal/constants.h"
#include "intrinsic/icon/reflexxes/internal/functional.h"
#include "intrinsic/icon/reflexxes/internal/math.h"
#include "intrinsic/icon/reflexxes/motion_polynomials.h"
#include "intrinsic/icon/reflexxes/outputs.h"
#include "intrinsic/icon/reflexxes/position_inputs.h"
#include "intrinsic/icon/reflexxes/position_outputs.h"
#include "intrinsic/icon/reflexxes/velocity_inputs.h"
#include "intrinsic/icon/reflexxes/velocity_outputs.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

namespace {
bool AlmostEqual(const double a, const double b) {
  return EpsilonEqual(a, b, kInputValueEpsilon);
}
}  // namespace

bool DoOutputsDeviateFromInputs(const Inputs& inputs, const Outputs& outputs) {
  if (!outputs.IsValidOutputAvailable()) {
    return true;
  }

  for (auto [dof_input, dof_output] :
       Zip(inputs.GetDOFs(), outputs.GetDOFs())) {
    if (dof_input.selected) {
      if (!(AlmostEqual(dof_output.new_position, dof_input.position) &&
            AlmostEqual(dof_output.new_velocity, dof_input.velocity) &&
            AlmostEqual(dof_output.new_acceleration, dof_input.acceleration))) {
        return true;
      }
    }
  }

  return false;
}

bool AreInputsCompatible(const Inputs& old_inputs, const Inputs& new_inputs) {
  for (auto [old_dof_input, new_dof_input] :
       Zip(old_inputs.GetDOFs(), new_inputs.GetDOFs())) {
    if (old_dof_input.selected != new_dof_input.selected) {
      return false;
    }

    // Skip if this isn't selected
    if (!old_dof_input.selected) {
      continue;
    }

    if (!AlmostEqual(old_dof_input.max_position, new_dof_input.max_position)) {
      return false;
    }

    if (!AlmostEqual(old_dof_input.min_position, new_dof_input.min_position)) {
      return false;
    }

    if (!AlmostEqual(old_dof_input.max_acceleration,
                     new_dof_input.max_acceleration)) {
      return false;
    }
    if (!AlmostEqual(old_dof_input.min_acceleration,
                     new_dof_input.min_acceleration)) {
      return false;
    }
    if (!AlmostEqual(old_dof_input.max_jerk, new_dof_input.max_jerk)) {
      return false;
    }
    if (!AlmostEqual(old_dof_input.min_jerk, new_dof_input.min_jerk)) {
      return false;
    }

    if (!AlmostEqual(old_dof_input.target_velocity,
                     new_dof_input.target_velocity)) {
      return false;
    }
  }

  // If we are passed in a min sync time, ensure that the old version and the
  // new version are the same but one cycle apart (as necessitated by one cycle
  // passing since the last call was made between the two inputs).
  if ((old_inputs.GetMinimumSynchronizationTime() > 0.0 ||
       new_inputs.GetMinimumSynchronizationTime() > 0.0) &&
      !AlmostEqual(old_inputs.GetMinimumSynchronizationTime() -
                       old_inputs.GetCycleTime(),
                   new_inputs.GetMinimumSynchronizationTime())) {
    return false;
  }

  return true;
}

bool AreInputsCompatible(const PositionInputs& old_inputs,
                         const PositionInputs& new_inputs) {
  if (!AreInputsCompatible(static_cast<const Inputs&>(old_inputs),
                           static_cast<const Inputs&>(new_inputs))) {
    return false;
  }

  for (auto [old_dof_input, new_dof_input] :
       Zip(old_inputs.GetDOFs(), new_inputs.GetDOFs())) {
    if (!AlmostEqual(old_dof_input.max_velocity, new_dof_input.max_velocity)) {
      return false;
    }
    if (!AlmostEqual(old_dof_input.min_velocity, new_dof_input.min_velocity)) {
      return false;
    }
    if (!AlmostEqual(old_dof_input.target_position,
                     new_dof_input.target_position)) {
      return false;
    }
  }

  return true;
}

bool ArePositionalLimitsBreached(const Inputs& inputs, const Outputs& outputs,
                                 double tolerance) {
  for (const auto [dof_input, dof_output] :
       Zip(inputs.GetDOFs(), outputs.GetDOFs())) {
    if (dof_output.max_position > (dof_input.max_position + tolerance) ||
        dof_output.min_position < (dof_input.min_position - tolerance)) {
      return true;
    }
  }
  return false;
}

bool AreVelocityLimitsBreached(const Inputs& inputs, const Outputs& outputs) {
  for (const auto [dof_input, dof_output] :
       Zip(inputs.GetDOFs(), outputs.GetDOFs())) {
    if (dof_output.max_velocity >
            dof_input.max_velocity * (1 + kValidSolutionEpsilon) ||
        dof_output.min_velocity <
            dof_input.min_velocity * (1 + kValidSolutionEpsilon)) {
      return true;
    }
  }
  return false;
}

bool AreAccelerationLimitsBreached(const Inputs& inputs,
                                   const Outputs& outputs) {
  for (const auto [dof_input, dof_output] :
       Zip(inputs.GetDOFs(), outputs.GetDOFs())) {
    if (dof_output.max_acceleration >
            dof_input.max_acceleration * (1 + kValidSolutionEpsilon) ||
        dof_output.min_acceleration <
            dof_input.min_acceleration * (1 + kValidSolutionEpsilon)) {
      return true;
    }
  }
  return false;
}

bool AreJerkLimitsBreached(const Inputs& inputs, const Outputs& outputs) {
  for (const auto [dof_input, dof_output] :
       Zip(inputs.GetDOFs(), outputs.GetDOFs())) {
    if (dof_output.max_jerk >
            dof_input.max_jerk * (1 + kValidSolutionEpsilon) ||
        dof_output.min_jerk <
            dof_input.min_jerk * (1 + kValidSolutionEpsilon)) {
      return true;
    }
  }
  return false;
}

namespace {
void ComputeExtremaData(Outputs& outputs) {
  for (Outputs::DOF& output_dof : outputs.GetDOFs()) {
    if (output_dof.motion_polynomials.GetSegmentCount() > 0) {
      // compute trajectory extrema until the last polynomial segment, ignoring
      // the last segment because it always contains an infinite time segment
      // after the goal is reached.
      MotionPolynomials::Extrema position_extrema =
          output_dof.motion_polynomials.GetPositionExtrema(
              output_dof.motion_polynomials.GetSegmentCount() - 1);
      output_dof.max_position = position_extrema.max_value;
      output_dof.min_position = position_extrema.min_value;
      output_dof.max_position_extrema_time = position_extrema.max_time;
      output_dof.min_position_extrema_time = position_extrema.min_time;

      MotionPolynomials::Extrema velocity_extrema =
          output_dof.motion_polynomials.GetVelocityExtrema();
      output_dof.max_velocity = velocity_extrema.max_value;
      output_dof.min_velocity = velocity_extrema.min_value;
      output_dof.max_velocity_extrema_time = velocity_extrema.max_time;
      output_dof.min_velocity_extrema_time = velocity_extrema.min_time;

      MotionPolynomials::Extrema acceleration_extrema =
          output_dof.motion_polynomials.GetAccelerationExtrema();
      output_dof.max_acceleration = acceleration_extrema.max_value;
      output_dof.min_acceleration = acceleration_extrema.min_value;
      output_dof.max_acceleration_extrema_time = acceleration_extrema.max_time;
      output_dof.min_acceleration_extrema_time = acceleration_extrema.min_time;

      MotionPolynomials::Extrema jerk_extrema =
          output_dof.motion_polynomials.GetJerkExtrema();
      output_dof.max_jerk = jerk_extrema.max_value;
      output_dof.min_jerk = jerk_extrema.min_value;
      output_dof.max_jerk_extrema_time = jerk_extrema.max_time;
      output_dof.min_jerk_extrema_time = jerk_extrema.min_time;
    }
  }
}
}  // namespace

void ComputeExtremaData(const PositionInputs& inputs,
                        PositionOutputs& outputs) {
  ComputeExtremaData(outputs);

  // Check if target position is exceeded.
  outputs.SetTrajectoryExceedsTargetPosition(false);
  for (auto [dof_input, dof_output] :
       Zip(inputs.GetDOFs(), outputs.GetDOFs())) {
    if (((dof_input.position < dof_input.target_position) &&
         (dof_output.max_position >
          (dof_input.target_position + kExceedingTargetPositionEpsilon))) ||
        ((dof_input.position > dof_input.target_position) &&
         (dof_output.min_position <
          (dof_input.target_position - kExceedingTargetPositionEpsilon)))) {
      outputs.SetTrajectoryExceedsTargetPosition(true);
    }
  }
}

void ComputeExtremaData(const VelocityInputs&, VelocityOutputs& outputs) {
  ComputeExtremaData(outputs);

  // compute position values at target velocity
  for (auto& dof_output : outputs.GetDOFs()) {
    // get the initial position of the last segment, since the last segment
    // extends to infinite time
    CHECK_GT(dof_output.motion_polynomials.GetSegmentCount(), 0);
    outputs.GetPositionValuesAtTargetVelocity()[dof_output.index] =
        dof_output.motion_polynomials.GetSegments().back().position;
  }
}

void CopySelectionVector(const Inputs& inputs, Outputs& outputs) {
  for (auto [dof_input, dof_output] :
       Zip(inputs.GetDOFs(), outputs.GetDOFs())) {
    dof_output.selected = dof_input.selected;
  }
}

void CopyNewStateToInputs(const Outputs& outputs, Inputs& inputs) {
  for (auto [dof_input, dof_output] :
       Zip(inputs.GetDOFs(), outputs.GetDOFs())) {
    dof_input.position = dof_output.new_position;
    dof_input.velocity = dof_output.new_velocity;
    dof_input.acceleration = dof_output.new_acceleration;
  }
}

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic
