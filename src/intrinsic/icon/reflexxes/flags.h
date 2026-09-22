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

#ifndef INTRINSIC_ICON_REFLEXXES_FLAGS_H_
#define INTRINSIC_ICON_REFLEXXES_FLAGS_H_

namespace intrinsic {
namespace reflexxes {

// Class containing flags to parameterize the execution of the Online Trajectory
// Generation algorithm.
// This class cannot be constructed directly, use VelocityFlags or PositionFlags
// instead.
struct Flags {
  // Enumeration whose values specify the synchronization method of the Online
  // Trajectory Generation algorithm.
  enum class SyncBehavior {
    // This is the default value. If it is possible to calculate a time-optimal
    // phase-synchronized (i.e., homothetic) trajectory, the algorithm will
    // generate it.
    kPhaseSynchronizationIfPossible = 0,
    // Even if it is possible to calculate a time-optimal phase-synchronized
    // trajectory, only a time-synchronized trajectory will be provided.
    kOnlyTimeSynchronization = 1,
    // Only phase-synchronized trajectories are allowed. If it is not possible
    // to calculate a phase-synchronized trajectory, an error will be returned
    // (but feasible, steady, and continuous output values will still be
    // computed).
    kOnlyPhaseSynchronization = 2,
    // When the boundary conditions are collinear, calculate a
    // phase-synchronized trajectory even if it is not time-optimal.
    kPhaseSynchronizationWhenCollinear = 3,
    // No synchronization will be performed, and all selected degrees of freedom
    // are treated independently.
    kNoSynchronization = 4,
  };

  // Enumeration whose values specify how positional limits are used:
  //
  // The trajectory generation algorithms can either prevent
  // - breaching positional limits or
  // - constraints for velocity, acceleration, jerk, etc.
  //
  // It is not possible to prevent from breaching both. In order to
  // guarantee executable motions, the trajectory generation algorithms of
  // the Reflexxes Motion Libraries do not change the constraints for
  // velocity, acceleration, and jerk.
  //
  // If positional limits will be breached by the currently computed
  // trajectory, the Reflexxes Motion Libraries provide three different
  // options to specify the behavior.
  //
  // -# Flags::kIgnore: ignore positional limits (default)
  // -# Flags::kErrorMsgOnly: return an error value but
  //           perform no action
  // -# Flags::kActivelyPrevent: return an error value
  //           and actively prevent from breaching the limits
  //
  // If an error value is returned because the currently computed trajectory
  // will exceed the positional limits, the user application will get
  // informed within one control cycle (typically one millisecond or less).
  // Depending on the application, the user can then implement a desired
  // behavior for such cases.
  enum class PositionalLimitsBehavior {
    // This is the default value. The values of Inputs::MaxPositionVector and
    // Inputs::MinPositionVector are ignored when calling
    // reflexxes::ComputePosition or reflexxes::ComputeVelocity
    kIgnore = 0,
    // If the trajectories will exceed the positional limits return an error
    // value, Status::kPositionalLimits.
    kErrorMsgOnly = 1,
    // Same behavior as in the case of kErrorMsgOnly, and in addition exceeding
    // the positional limits will be actively prevented. This active
    // prevention will also be performed in case the desired state of motion has
    // already been reached.
    kActivelyPrevent = 2,
  };

  // Behavior when the constraints are invalid
  enum class InvalidConstraintsBehavior {
    // Deselect the degree of freedom and report no error.
    kDeselectDofWithoutErrorMsg = 0,
    // Report error and quit online trajectory generation.
    kReportErrorAndQuit = 1,
  };

  // Behavior when the scale of input values is too large.
  enum class InvalidScaleOfInputValuesBehavior {
    // Ignore.
    kIgnore = 0,
    // Ignore but report error.
    kErrorMsgOnly = 1,
    // Report error and quit online trajectory generation.
    kReportErrorAndQuit = 2,
  };

  // Checks if the object values are equal to that of the input object.
  bool operator==(const Flags& flags) const {
    return (
        (synchronization_behavior == flags.synchronization_behavior) &&
        (positional_limits_behavior == flags.positional_limits_behavior) &&
        (invalid_constraints_behavior == flags.invalid_constraints_behavior) &&
        (invalid_scale_of_input_values_behavior ==
         flags.invalid_scale_of_input_values_behavior));
  }

  // Checks if the object values are unequal to that of the input object.
  bool operator!=(const Flags& flags) const { return !(*this == flags); }

  SyncBehavior synchronization_behavior;
  PositionalLimitsBehavior positional_limits_behavior;
  InvalidConstraintsBehavior invalid_constraints_behavior;
  InvalidScaleOfInputValuesBehavior invalid_scale_of_input_values_behavior;

 protected:
  // Since the child classes set these differently, force setting all members in
  // the constructor.
  Flags(
      SyncBehavior synchronization_behavior,
      PositionalLimitsBehavior positional_limits_behavior,
      InvalidConstraintsBehavior invalid_constraints_behavior,
      InvalidScaleOfInputValuesBehavior invalid_scale_of_input_values_behavior)
      : synchronization_behavior(synchronization_behavior),
        positional_limits_behavior(positional_limits_behavior),
        invalid_constraints_behavior(invalid_constraints_behavior),
        invalid_scale_of_input_values_behavior(
            invalid_scale_of_input_values_behavior) {}
};

}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_FLAGS_H_
