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

#include "intrinsic/icon/reflexxes/reflexxes_api.h"

#include <algorithm>
#include <cmath>
#include <optional>

#include "absl/algorithm/container.h"
#include "absl/log/check.h"
#include "intrinsic/icon/reflexxes/flags.h"
#include "intrinsic/icon/reflexxes/inputs.h"
#include "intrinsic/icon/reflexxes/internal/calc_min_execution_time_utils.h"
#include "intrinsic/icon/reflexxes/internal/constants.h"
#include "intrinsic/icon/reflexxes/internal/functional.h"
#include "intrinsic/icon/reflexxes/internal/input_output_utility_functions.h"
#include "intrinsic/icon/reflexxes/internal/math.h"
#include "intrinsic/icon/reflexxes/internal/position.h"
#include "intrinsic/icon/reflexxes/internal/velocity.h"
#include "intrinsic/icon/reflexxes/motion_polynomials.h"
#include "intrinsic/icon/reflexxes/outputs.h"
#include "intrinsic/icon/reflexxes/position_flags.h"
#include "intrinsic/icon/reflexxes/position_inputs.h"
#include "intrinsic/icon/reflexxes/position_outputs.h"
#include "intrinsic/icon/reflexxes/status.h"
#include "intrinsic/icon/reflexxes/velocity_flags.h"
#include "intrinsic/icon/reflexxes/velocity_inputs.h"
#include "intrinsic/icon/reflexxes/velocity_outputs.h"

namespace intrinsic {
namespace reflexxes {

State::State(const int num_dofs, const double cycle_time)
    : position_state(num_dofs, cycle_time),
      velocity_state(num_dofs, cycle_time),
      fallback_velocity_state(num_dofs, cycle_time) {}

int State::GetNumberOfDOFs() const {
  return position_state.inputs.GetNumberOfDOFs();
}

double State::GetCycleTime() const {
  return position_state.inputs.GetCycleTime();
}

namespace internal {
namespace {

// Gets the next state of motion based on the current state and the new inputs.
// Prototype this here as some of the fallback functions use it.
template <typename InputT, typename FlagsT, typename OutputT, typename StateT>
Status GetNextStateOfMotion(const InputT& inputs, const FlagsT& flags,
                            OutputT& outputs, StateT& state);

bool CheckForFiniteness(const Inputs& input) {
  for (const Inputs::DOF& dof : input.GetDOFs()) {
    if (!isfinite(dof.position) || !isfinite(dof.velocity) ||
        !isfinite(dof.acceleration) || !isfinite(dof.target_position) ||
        !isfinite(dof.target_velocity) || isnan(dof.max_position) ||
        !isfinite(dof.max_velocity) || !isfinite(dof.max_acceleration) ||
        !isfinite(dof.max_jerk) || isnan(dof.min_position) ||
        !isfinite(dof.min_velocity) || !isfinite(dof.min_acceleration) ||
        !isfinite(dof.min_jerk)) {
      return false;
    }
  }

  return isfinite(input.GetMinimumSynchronizationTime());
}

// Executes the fallback strategy for ComputePosition.
// Tries to meet target velocity in the next cycle with ComputeVelocity.
void FallbackStrategy(const PositionInputs& inputs, const PositionFlags& flags,
                      PositionOutputs& outputs,
                      internal::VelocityState& fallback_velocity_state) {
  VelocityInputs velocity_inputs(inputs);
  VelocityFlags velocity_flags(flags);
  VelocityOutputs velocity_output(velocity_inputs.GetNumberOfDOFs(),
                                  velocity_inputs.GetCycleTime());

  for (auto [input_dof, velocity_dof] :
       Zip(inputs.GetDOFs(), velocity_inputs.GetDOFs())) {
    velocity_dof.target_velocity =
        flags.keep_current_velocity_in_case_of_fallback_strategy
            ? input_dof.target_velocity
            : input_dof.alt_target_velocity;
  }

  // If phase synchronization is not required, then choose no synchronization.
  if (flags.synchronization_behavior !=
      Flags::SyncBehavior::kOnlyPhaseSynchronization) {
    velocity_flags.synchronization_behavior =
        Flags::SyncBehavior::kNoSynchronization;
  }

  GetNextStateOfMotion(velocity_inputs, velocity_flags, velocity_output,
                       fallback_velocity_state);

  outputs = velocity_output;
}

// Executes the fallback strategy for ComputePosition.
// Tries to meet target velocity in the next cycle by continuing with current
// acceleration through the next cycle.
void FallbackStrategy(const VelocityInputs& inputs, const VelocityFlags&,
                      VelocityOutputs& outputs, internal::VelocityState&) {
  for (auto [dof_input, dof_output] :
       Zip(inputs.GetDOFs(), outputs.GetDOFs())) {
    dof_output.new_position =
        dof_input.position + dof_input.velocity * inputs.GetCycleTime() +
        0.5 * dof_input.acceleration * internal::Power2(inputs.GetCycleTime());
    dof_output.new_velocity =
        dof_input.velocity + dof_input.acceleration * inputs.GetCycleTime();
    dof_output.new_acceleration = dof_input.acceleration;
  }
}

// Generates a velocity strategy to bring all DOFs to a halt.
void AsynchronousStopMotion(const Inputs& inputs, const Flags& flags,
                            Outputs& outputs,
                            internal::VelocityState& fallback_velocity_state) {
  VelocityInputs velocity_input(inputs);
  for (Inputs::DOF& dof : velocity_input.GetDOFs()) {
    dof.target_velocity = 0.0;
  }

  VelocityFlags velocity_flags(flags);
  velocity_flags.synchronization_behavior =
      Flags::SyncBehavior::kNoSynchronization;

  VelocityOutputs velocity_output(velocity_input.GetNumberOfDOFs(),
                                  velocity_input.GetCycleTime());

  internal::GetNextStateOfMotion(velocity_input, velocity_flags,
                                 velocity_output, fallback_velocity_state);

  outputs = velocity_output;
}

// Returns whether the position at target velocity breaches the positional
// limits.
bool DoesPositionAtTargetVelocityBreachPositionalLimits(
    const Inputs& inputs, const Outputs& outputs,
    double tolerance = kPositionLimitBreachedThreshold) {
  for (auto [dof_input, dof_output] :
       Zip(inputs.GetDOFs(), outputs.GetDOFs())) {
    // Get the final position by getting the last segment before the capping
    // "infinite segment" and finding its end time.
    double end_time = 0.0;
    const auto& segments = dof_output.motion_polynomials.GetSegments();
    auto iter = std::find_if(segments.rbegin(), segments.rend(),
                             [](const auto& segment) {
                               return segment.end_time != internal::kInfinity;
                             });

    if (iter != segments.rend()) {
      end_time = iter->end_time;
    }

    MotionPolynomials::MotionState end_motion_state =
        dof_output.motion_polynomials.GetStateOfMotionAtTime(end_time);

    // Check for breach.
    if (end_motion_state.position > dof_input.max_position + tolerance ||
        end_motion_state.position < dof_input.min_position - tolerance) {
      return true;
    }
  }
  return false;
}

bool DoesStopMotionFromNewStateBreachPositionalLimits(const Inputs& inputs,
                                                      const Flags& flags,
                                                      const Outputs& outputs) {
  // Copy new state to current state and set a target velocity of 0.
  VelocityInputs velocity_inputs(inputs);
  internal::CopyNewStateToInputs(outputs, velocity_inputs);
  for (Inputs::DOF& dof : velocity_inputs.GetDOFs()) {
    dof.target_velocity = 0.0;
  }

  VelocityFlags velocity_flags = flags;
  velocity_flags.synchronization_behavior =
      Flags::SyncBehavior::kNoSynchronization;

  // Reset
  VelocityOutputs velocity_outputs(velocity_inputs.GetNumberOfDOFs(),
                                   velocity_inputs.GetCycleTime());
  internal::VelocityState velocity_state(velocity_inputs.GetNumberOfDOFs(),
                                         velocity_inputs.GetCycleTime());

  // Generate stop motion and check if it breaches.
  internal::GetNextStateOfMotion(velocity_inputs, velocity_flags,
                                 velocity_outputs, velocity_state);

  return DoesPositionAtTargetVelocityBreachPositionalLimits(velocity_inputs,
                                                            velocity_outputs);
}

template <typename InputT, typename FlagsT, typename OutputT, typename StateT>
std::optional<Status> CheckInputs(const InputT& inputs, const FlagsT& flags,
                                  const OutputT& outputs, const StateT& state) {
  if (inputs.GetNumberOfDOFs() != outputs.GetNumberOfDOFs()) {
    return Status::kErrorNumberOfDofs;
  }

  if (inputs.GetNumberOfDOFs() != state.inputs.GetNumberOfDOFs()) {
    return Status::kErrorNumberOfDofs;
  }

  if (inputs.GetCycleTime() != outputs.GetCycleTime()) {
    return Status::kErrorCycleTime;
  }

  if (!CheckForFiniteness(inputs)) {
    return Status::kErrorInvalidInputValues;
  }

  if (!inputs.CheckValidityOfConstraints() &&
      flags.invalid_constraints_behavior ==
          Flags::InvalidConstraintsBehavior::kReportErrorAndQuit) {
    return Status::kErrorInvalidConstraints;
  }

  return std::nullopt;
}

// Scale the outputs with the inverse of the input scaling.
void ScaleOutput(VelocityOutputs& outputs,
                 const MaxDOFFixedVector<double>& scaling_factors) {
  for (VelocityOutputs::DOF& dof : outputs.GetDOFs()) {
    dof.Scale(1.0 / scaling_factors[dof.index]);
    outputs.GetPositionValuesAtTargetVelocity()[dof.index] /=
        scaling_factors[dof.index];
  }
}

void ScaleOutput(PositionOutputs& outputs,
                 const MaxDOFFixedVector<double>& scaling_factors) {
  for (PositionOutputs::DOF& dof : outputs.GetDOFs()) {
    dof.Scale(1.0 / scaling_factors[dof.index]);
  }
}

// Runs Step1 & Step2 of the desired algorithm and processes the results.
template <typename InputT, typename FlagsT, typename OutputT>
std::optional<Status> RunSteps(const FlagsT& flags, InputT& inputs,
                               OutputT& outputs, bool force_ignore_phase_sync) {
  InputT original_input_backup = inputs;

  MaxDOFFixedVector<double> scaling_factors = inputs.ScaleIfNecessary();

  // Based on the scaling, figure out if any DOFs need to be deselected.
  inputs.DeselectInvalidDofs();
  CopySelectionVector(inputs, outputs);

  if (!CalcMinExecutionTime(inputs, flags, outputs)) {
    // Note the check after the fact, and not before.  Despite failing a
    // validity check, if step 1 is valid the algorithm proceeds.
    if (!inputs.IsValid()) {
      return Status::kErrorInvalidInputValues;
    }
    return Status::kErrorExecutionTimeCalculation;
  }

  // Make sure to force disabling phase sync.
  if (force_ignore_phase_sync) {
    outputs.DisablePhaseSync();
  }

  if (outputs.GetSyncTime() > kMaxAllowedSyncTimeSeconds) {
    return Status::kErrorExecutionTimeTooBig;
  }
  // check if only phase synchronization is requested
  // if yes, check if it is possible
  if ((flags.synchronization_behavior ==
       Flags::SyncBehavior::kOnlyPhaseSynchronization) &&
      !outputs.IsPhaseSyncEnabled()) {
    return Status::kErrorNoPhaseSynchronization;
  }

  // A precondition for step 2 is scaling the limits in case of phase
  // synchronization.
  if (outputs.IsPhaseSyncEnabled() && !force_ignore_phase_sync) {
    ScaleLimitsForPhaseSync(inputs, outputs);
  }

  // Reflexxes step 2.
  // Perform a phase/time synchronization of the degrees of freedom.
  const bool was_dofs_sync_successful = SynchronizeDOFs(inputs, flags, outputs);

  // We only run the next checks if `SynchronizeDOFs` succeeded.
  bool are_positional_limits_breached = true;
  if (was_dofs_sync_successful) {
    // Compute position extrema now that new trajectories are available.
    ComputeExtremaData(inputs, outputs);

    // If the user did not request to ignore positional limits check for a
    // breach. It's important to do this *before* scaling the outputs back up,
    // because we compare against the scaled-down limits!
    are_positional_limits_breached =
        (flags.positional_limits_behavior !=
             Flags::PositionalLimitsBehavior::kIgnore &&
         ArePositionalLimitsBreached(inputs, outputs));
  }

  if (!was_dofs_sync_successful || are_positional_limits_breached) {
    // If SynchronizeDOFs failed and current call does not enforce ignoring
    // phase sync, try to recover by falling back to time synchronization.
    // Recursively calls `RunSteps` and forces to ignore phase synchronization
    // in the next attempt.
    // Reflexxes identifies collinear inputs and attempts to build phase
    // synchronized trajectories with them. However, in the current version the
    // computed phase synchronized trajectories could violate positional limits.
    // Thus, here we attempt to build only a time synchronized solution, in a
    // similar fashion as when the phase synchronization was not successful.
    if (!force_ignore_phase_sync) {
      inputs = original_input_backup;
      return RunSteps(flags, inputs, outputs, /*force_ignore_phase_sync=*/true);
    }
    return (!was_dofs_sync_successful ? Status::kErrorSynchronization
                                      : Status::kErrorPositionalLimits);
  }

  // check if only phase synchronization is requested
  // if yes, check if it is still possible after step2
  if ((flags.synchronization_behavior ==
       Flags::SyncBehavior::kOnlyPhaseSynchronization) &&
      !outputs.IsPhaseSyncEnabled()) {
    return Status::kErrorNoPhaseSynchronization;
  }

  // Due to inconsistencies in phase synchronization, we need to examine the
  // velocity and acceleration trajectories for limit breaches. In case of
  // velocity/acceleration limit breaches, we fall back to standard
  // time-synchronization, however only if 'kPhaseSynchronizationWhenCollinear'
  // is selected. For more details, see b/223128905.
  if (outputs.IsPhaseSyncEnabled() &&
      flags.synchronization_behavior ==
          Flags::SyncBehavior::kPhaseSynchronizationWhenCollinear) {
    if (AreVelocityLimitsBreached(original_input_backup, outputs) ||
        AreAccelerationLimitsBreached(original_input_backup, outputs) ||
        AreJerkLimitsBreached(original_input_backup, outputs)) {
      inputs = original_input_backup;
      return RunSteps(flags, inputs, outputs, /*force_ignore_phase_sync=*/true);
    }
  }

  ScaleOutput(outputs, scaling_factors);

  return std::nullopt;
}

bool IsReplanningNecessary(const VelocityInputs& inputs,
                           const VelocityFlags& flags,
                           const VelocityState& state) {
  return state.flags != flags || !AreInputsCompatible(state.inputs, inputs) ||
         DoOutputsDeviateFromInputs(inputs, state.outputs);
}

bool IsReplanningNecessary(const PositionInputs& inputs,
                           const PositionFlags& flags,
                           const PositionState& state) {
  return state.flags != flags || !AreInputsCompatible(state.inputs, inputs) ||
         DoOutputsDeviateFromInputs(inputs, state.outputs) ||
         (state.outputs.GetStatus() == Status::kFinalStateReached &&
          flags.final_motion_behavior ==
              PositionFlags::FinalMotionBehavior::kRecomputeTrajectory);
}

bool AreTargetsValid(const PositionInputs& inputs, bool check_position) {
  auto pos_valid = [check_position](const Inputs::DOF& dof) {
    return !check_position ||
           IsInRange({dof.min_position, dof.max_position}, dof.target_position);
  };
  auto vel_valid = [](const Inputs::DOF& dof) {
    return IsInRange({dof.min_velocity, dof.max_velocity}, dof.target_velocity);
  };

  return absl::c_all_of(inputs.GetDOFs(), [&](const Inputs::DOF& dof) {
    return !dof.selected || (pos_valid(dof) && vel_valid(dof));
  });
}

bool AreTargetsValid(const VelocityInputs&, bool) { return true; }

// Gets the next state of motion based on the current state and the new inputs.
template <typename InputT, typename FlagsT, typename OutputT, typename StateT>
Status GetNextStateOfMotion(const InputT& inputs, const FlagsT& flags,
                            OutputT& outputs, StateT& state) {
  // Check whether target state is valid.
  if (!AreTargetsValid(inputs,
                       flags.positional_limits_behavior ==
                           Flags::PositionalLimitsBehavior::kActivelyPrevent)) {
    outputs.SetStatus(Status::kErrorInvalidTargetState);
    return outputs.GetStatus();
  }

  // Check whether the scale of input values is valid.
  if (!inputs.IsScaleOfInputsValid() &&
      flags.invalid_scale_of_input_values_behavior ==
          Flags::InvalidScaleOfInputValuesBehavior::kReportErrorAndQuit) {
    outputs.SetStatus(Status::kErrorInvalidScaleOfInputValues);
    return outputs.GetStatus();
  }

  // Check whether any planning/re-planning is necessary.
  std::optional<Status> opt_status;
  if (IsReplanningNecessary(inputs, flags, state)) {
    // Use the state.inputs as a scratch space to scale.
    state.inputs = inputs;

    // Reset the outputs
    state.outputs = OutputT(inputs.GetNumberOfDOFs(), inputs.GetCycleTime());

    opt_status = RunSteps(flags, state.inputs, state.outputs,
                          /*force_ignore_phase_sync=*/false);
    state.inputs = inputs;
    state.flags = flags;
  } else {
    // Check for existing positional limits error and use the same status
    // message if it's there.
    if (state.outputs.GetStatus() == Status::kErrorPositionalLimits) {
      opt_status = Status::kErrorPositionalLimits;
    }
    // Check if the minimum synchronization time needs to be updated.
    if (state.inputs.GetMinimumSynchronizationTime() > 0.0) {
      state.inputs.SetMinimumSynchronizationTime(
          state.inputs.GetMinimumSynchronizationTime() -
          state.inputs.GetCycleTime());
    }
  }

  if (flags.invalid_scale_of_input_values_behavior ==
          Flags::InvalidScaleOfInputValuesBehavior::kErrorMsgOnly &&
      !inputs.IsScaleOfInputsValid()) {
    opt_status = Status::kErrorInvalidScaleOfInputValues;
  }

  // Compute next state if there is no error that prevents it.
  if (!opt_status || *opt_status == Status::kErrorPositionalLimits ||
      *opt_status == Status::kErrorInvalidScaleOfInputValues) {
    Status output_status = state.outputs.ComputeNextStateOfMotion();
    // If there was no error message then use output status message.
    if (!opt_status) {
      opt_status = output_status;
    }
  }

  state.outputs.SetStatus(*opt_status);
  outputs = state.outputs;
  return *opt_status;
}

template <typename InputT, typename FlagsT, typename OutputT, typename StateT>
Status Compute(const InputT& inputs, const FlagsT& flags, OutputT& outputs,
               StateT& state,
               internal::VelocityState& fallback_velocity_state) {
  std::optional<Status> opt_status = CheckInputs(inputs, flags, outputs, state);

  if (opt_status) {
    outputs.SetStatus(*opt_status);
    return *opt_status;
  }

  Status status = GetNextStateOfMotion(inputs, flags, outputs, state);

  if (status != Status::kWorking && status != Status::kFinalStateReached) {
    if (status == Status::kErrorPositionalLimits) {
      if (flags.positional_limits_behavior ==
          Flags::PositionalLimitsBehavior::kActivelyPrevent) {
        // User has requested that breaches be prevented, so use a fallback
        // strategy to stop the motion.
        outputs = OutputT(inputs.GetNumberOfDOFs(), inputs.GetCycleTime());
        AsynchronousStopMotion(inputs, flags, outputs, fallback_velocity_state);
      }
    } else {
      // Fallback unless the scale is invalid or the execution time is too big.
      if (status != Status::kErrorInvalidScaleOfInputValues &&
          status != Status::kErrorExecutionTimeTooBig) {
        FallbackStrategy(inputs, flags, outputs, fallback_velocity_state);
      }
    }
  }

  // Check if stop motion from the new state breaches positional limits, and
  // then stop it if the user has requested it.
  if (flags.positional_limits_behavior !=
      Flags::PositionalLimitsBehavior::kIgnore) {
    if (DoesStopMotionFromNewStateBreachPositionalLimits(inputs, flags,
                                                         outputs)) {
      status = Status::kErrorPositionalLimits;
      if (flags.positional_limits_behavior ==
          Flags::PositionalLimitsBehavior::kActivelyPrevent) {
        outputs = OutputT(inputs.GetNumberOfDOFs(), inputs.GetCycleTime());
        AsynchronousStopMotion(inputs, flags, outputs, fallback_velocity_state);
      }
    }
  }

  outputs.SetStatus(status);
  return status;
}

}  // namespace
}  // namespace internal

Status ComputePosition(const PositionInputs& inputs, const PositionFlags& flags,
                       PositionOutputs& outputs, State& state) {
  return internal::Compute(inputs, flags, outputs, state.position_state,
                           state.fallback_velocity_state);
}

Status ComputeVelocity(const VelocityInputs& inputs, const VelocityFlags& flags,
                       VelocityOutputs& outputs, State& state) {
  return internal::Compute(inputs, flags, outputs, state.velocity_state,
                           state.fallback_velocity_state);
}

bool DoesStopMotionFromNewStateBreachPositionalLimits(
    const VelocityInputs& inputs, const VelocityFlags& flags,
    const VelocityOutputs& outputs) {
  CHECK_EQ(inputs.GetNumberOfDOFs(), outputs.GetNumberOfDOFs());
  CHECK_EQ(inputs.GetCycleTime(), outputs.GetCycleTime());
  CHECK_EQ(internal::CheckForFiniteness(inputs), true);
  if (flags.invalid_constraints_behavior ==
      Flags::InvalidConstraintsBehavior::kReportErrorAndQuit) {
    CHECK_EQ(inputs.CheckValidityOfConstraints(), true);
  }

  return internal::DoesStopMotionFromNewStateBreachPositionalLimits(
      inputs, flags, outputs);
}

}  // namespace reflexxes
}  // namespace intrinsic
