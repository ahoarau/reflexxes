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

#include "intrinsic/icon/reflexxes/internal/synchronize_dofs_base.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <utility>

#include "absl/functional/function_ref.h"
#include "absl/log/check.h"
#include "intrinsic/icon/reflexxes/flags.h"
#include "intrinsic/icon/reflexxes/inputs.h"
#include "intrinsic/icon/reflexxes/internal/constants.h"
#include "intrinsic/icon/reflexxes/internal/decision_tree_utility_functions.h"
#include "intrinsic/icon/reflexxes/internal/functional.h"
#include "intrinsic/icon/reflexxes/internal/math.h"
#include "intrinsic/icon/reflexxes/internal/motion_state.h"
#include "intrinsic/icon/reflexxes/motion_polynomials.h"
#include "intrinsic/icon/reflexxes/outputs.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {
namespace {

// Sets the output trajectory polynomial to be an extrapolation
// of the current state for the unselected dof
void SetExtrapolatedTrajectoryWithCurrentState(const Inputs::DOF& dof_input,
                                               Outputs::DOF& dof_output) {
  // reset motion polynomials
  dof_output.motion_polynomials = MotionPolynomials();
  dof_output.motion_polynomials.AddSegment(
      {.position = dof_input.position,
       .velocity = dof_input.velocity,
       .acceleration = dof_input.acceleration,
       .jerk = 0.0,
       .start_time = 0.0,
       .end_time = kInfinity});
}

// Tries to compute a timescale from the given extrema and potential min/max
// values, but only if the max value exceeds the max extrema, or the min is less
// than the min extrema, and only if the corresponding "check" flags are set.
// If a value is found, it is returned, if not -1.0 is returned.
double ComputeTimeScale(const double phase_sync_scale, const Range& extrema,
                        const Range& min_max, bool check_min, bool check_max,
                        const std::function<double(double)>& ratio_transform) {
  constexpr auto apply_epsilon = [](double original) {
    return original * (1.0 + kRelPhaseSyncEpsilon) - kAbsPhaseSyncEpsilon;
  };

  double time_scale = -1.0;
  Range scaled_extrema = extrema;
  scaled_extrema.min *= phase_sync_scale;
  scaled_extrema.max *= phase_sync_scale;
  if (phase_sync_scale <= 0) {
    std::swap(scaled_extrema.min, scaled_extrema.max);
  }

  if (check_max && scaled_extrema.max > apply_epsilon(min_max.max)) {
    time_scale =
        std::max(time_scale, ratio_transform(scaled_extrema.max / min_max.max));
  }

  if (check_min && scaled_extrema.min < apply_epsilon(min_max.min)) {
    time_scale =
        std::max(time_scale, ratio_transform(scaled_extrema.min / min_max.min));
  }

  return time_scale;
}

double ComputeJerkBasedTimeScale(const double phase_sync_scale,
                                 const Range& extrema, const Range& min_max,
                                 bool check_min, bool check_max) {
  return ComputeTimeScale(phase_sync_scale, extrema, min_max, check_min,
                          check_max,
                          static_cast<double (*)(double)>(&std::cbrt));
}

double ComputeAccelerationBasedTimeScale(const double phase_sync_scale,
                                         const Range& extrema,
                                         const Range& min_max) {
  return ComputeTimeScale(phase_sync_scale, extrema, min_max, true, true,
                          static_cast<double (*)(double)>(&std::sqrt));
}

double ComputeVelocityBasedTimeScale(const double phase_sync_scale,
                                     const Range& extrema,
                                     const Range& min_max) {
  return ComputeTimeScale(phase_sync_scale, extrema, min_max, true, true,
                          [](double x) { return x; });
}

// Calculates the timescale from the given inputs, the phase sync scale, and the
// calculated extrema. Returns the timescale if found, or -1.0 otherwise.
double CalcTimeScaleFromLimits(
    const Inputs::DOF& dof_input, const Inputs::DOF& dof_input_phase_sync,
    const double phase_sync_scale, const bool check_velocity_limits,
    const MotionPolynomials::Extrema& vel_extrema_phase_sync,
    const MotionPolynomials::Extrema& acc_extrema_phase_sync) {
  double time_scale = -1.;

  time_scale = std::max(
      time_scale,
      ComputeAccelerationBasedTimeScale(
          phase_sync_scale,
          {acc_extrema_phase_sync.min_value, acc_extrema_phase_sync.max_value},
          {dof_input.min_acceleration, dof_input.max_acceleration}));

  if (check_velocity_limits) {
    time_scale = std::max(
        time_scale, ComputeVelocityBasedTimeScale(
                        phase_sync_scale,
                        {vel_extrema_phase_sync.min_value,
                         vel_extrema_phase_sync.max_value},
                        {dof_input.min_velocity, dof_input.max_velocity}));
  }

  const bool check_max_jerk =
      (acc_extrema_phase_sync.max_value > dof_input.acceleration) ||
      (acc_extrema_phase_sync.min_value < 0.0);
  const bool check_min_jerk =
      (acc_extrema_phase_sync.min_value < dof_input_phase_sync.acceleration) ||
      (acc_extrema_phase_sync.max_value > 0.0);

  time_scale = std::max(
      time_scale,
      ComputeJerkBasedTimeScale(
          phase_sync_scale,
          {dof_input_phase_sync.min_jerk, dof_input_phase_sync.max_jerk},
          {dof_input.min_jerk, dof_input.max_jerk}, check_min_jerk,
          check_max_jerk));

  return time_scale;
}

void PhaseSyncScale(const Inputs& inputs, const bool phase_sync_when_collinear,
                    const bool check_velocity_limits, Outputs& outputs) {
  if (!outputs.IsPhaseSyncEnabled()) {
    return;
  }

  const int phase_sync_index = outputs.GetPhaseSyncDOFIndex();
  const auto& dof_input_phase_sync = inputs.GetDOFs()[phase_sync_index];
  auto& dof_output_phase_sync = outputs.GetDOFs()[phase_sync_index];

  // Retrieve MotionPolynomials and its extrema for phase sync dof
  MotionPolynomials motion_poly_phase_sync =
      dof_output_phase_sync.motion_polynomials;
  MotionPolynomials::Extrema vel_extrema_phase_sync =
      motion_poly_phase_sync.GetVelocityExtrema();
  MotionPolynomials::Extrema acc_extrema_phase_sync =
      motion_poly_phase_sync.GetAccelerationExtrema();
  // time scale factor for phase synchronized trajectories when collinear
  double time_scale = 1.0;
  // Check whether the phase sync scaling factors are valid
  for (auto [dof_input, dof_output] :
       Zip(inputs.GetDOFs(), outputs.GetDOFs())) {
    // for all selected dofs except the phase sync dof
    if (dof_input.selected && (dof_input.index != phase_sync_index)) {
      double dof_time_scale = CalcTimeScaleFromLimits(
          dof_input, dof_input_phase_sync, dof_output.phase_sync_scale,
          check_velocity_limits, vel_extrema_phase_sync,
          acc_extrema_phase_sync);
      if (dof_time_scale >= 0. && !phase_sync_when_collinear) {
        outputs.DisablePhaseSync();
        return;
      }

      time_scale = std::max(time_scale, dof_time_scale);
    }
  }

  // if phase sync when collinear flag is on, time scale the motion
  // polynomials
  if (phase_sync_when_collinear) {
    motion_poly_phase_sync.TimeScale(time_scale);
    outputs.SetSyncTime(outputs.GetSyncTime() * time_scale);
  }

  // loop through all dofs and overwrite their motion polynomials
  for (auto [dof_input, dof_output] :
       Zip(inputs.GetDOFs(), outputs.GetDOFs())) {
    Outputs::DOF::SubStep& substep_output = dof_output.step2;
    if (dof_input.selected) {
      // update the profile error element
      substep_output.result = true;

      // set applied profile element
      substep_output.applied_profile =
          dof_output_phase_sync.step2.applied_profile;
      // set applied profile trace from the phase sync dof
      if (dof_output.index != phase_sync_index) {
        substep_output.applied_profile_trace =
            dof_output_phase_sync.step2.applied_profile_trace;
      }
      // scale output polynomials
      dof_output.motion_polynomials = motion_poly_phase_sync;
      dof_output.motion_polynomials.ScaleWithNewInitialValues(
          dof_output.phase_sync_scale, dof_input.position, dof_input.velocity,
          dof_input.acceleration);

      // enforce the last segment to have zero acceleration and
      // target position and target velocity as coefficients
      CHECK_GT(dof_output.motion_polynomials.GetSegmentCount(), 0);
      MotionPolynomials::Segment& last_segment =
          dof_output.motion_polynomials.GetSegments().back();
      double ptrgt = check_velocity_limits ? dof_input.target_position
                                           : last_segment.position;

      last_segment.position = ptrgt;
      last_segment.velocity = dof_input.target_velocity;
      last_segment.acceleration = last_segment.jerk = 0.0;
    } else {
      // dof not selected
      // return success
      // update the profile error element
      substep_output.result = true;
      // extrapolate with current state
      SetExtrapolatedTrajectoryWithCurrentState(dof_input, dof_output);
    }
  }
}

}  // namespace

bool SynchronizeDOFs(
    const Inputs& inputs, const SyncDOFsFlags& flags,
    const absl::FunctionRef<SyncDOFsMotionState(const SyncDOFsMotionState&)>
        exec_decision_tree,
    Outputs& outputs) {
  for (auto [dof_input, dof_output] :
       Zip(inputs.GetDOFs(), outputs.GetDOFs())) {
    // Run Step2 for all dofs irrespective of whether phase sync is possible or
    // not

    if (!dof_input.selected) {
      // If the DOF is not selected, skip it and continue until the other
      // selected DOFs succeed or fail.
      dof_output.step2.result = true;

      // extrapolate with current state
      SetExtrapolatedTrajectoryWithCurrentState(dof_input, dof_output);
      continue;
    }

    StateBase state_base(
        dof_input, dof_output, dof_output.step2, &dof_output.motion_polynomials,
        dof_output.execution_time, flags.get_into_boundaries_fast);
    SyncDOFsMotionState s =
        exec_decision_tree(InitializeMotionState(state_base, dof_input));

    if (IsFailure(s)) {
      // an error occurred so quit deciding for the rest of the dofs
      // and notify failure
      return false;
    }

    WriteStateToSubStepOutput(s, dof_output.step2);
  }

  // check if phase synchronization is to be done
  PhaseSyncScale(inputs,
                 flags.sync_behavior ==
                     Flags::SyncBehavior::kPhaseSynchronizationWhenCollinear,
                 flags.use_positional_checks, outputs);

  return true;
}

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic
