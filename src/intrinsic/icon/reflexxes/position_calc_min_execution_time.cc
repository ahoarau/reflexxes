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

#include <algorithm>

#include "absl/algorithm/container.h"
#include "absl/log/check.h"
#include "intrinsic/icon/reflexxes/constants_external.h"
#include "intrinsic/icon/reflexxes/flags.h"
#include "intrinsic/icon/reflexxes/inputs.h"
#include "intrinsic/icon/reflexxes/internal/calc_min_execution_time_utils.h"
#include "intrinsic/icon/reflexxes/internal/constants.h"
#include "intrinsic/icon/reflexxes/internal/position.h"
#include "intrinsic/icon/reflexxes/outputs.h"
#include "intrinsic/icon/reflexxes/position_flags.h"
#include "intrinsic/icon/reflexxes/position_inputs.h"
#include "intrinsic/icon/reflexxes/position_outputs.h"
#include "intrinsic/icon/reflexxes/profile.h"
#include "intrinsic/util/fixed_vector.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

namespace {

// Add an "override" to the profile trace, resulting from the chosen profile
// being overridden
inline void OverrideProfileTrace(Outputs::DOF::SubStep& substep) {
  if constexpr (kEnableTracing) {
    CHECK(substep.applied_profile_trace_size < kTraceSize);
    substep.applied_profile_trace[substep.applied_profile_trace_size++] =
        Profile::kOverridden;
  }
}

// Ensure inoperative time range is at least the minimum execution time, and
// that the end inoperative time is larger than the beginning time.
void FixExecutionTimes(PositionOutputs& output_params) {
  for (auto& dof_output : output_params.GetDOFs()) {
    if (!dof_output.selected) {
      continue;
    }

    dof_output.min_execution_time += kStep1TimeEpsilon;

    // check if phase sync is on,
    // if yes ignore the rest
    if (output_params.IsPhaseSyncEnabled()) {
      continue;
    }

    if ((dof_output.inoperative_begin_execution_time != kInfinity) &&
        (dof_output.inoperative_end_execution_time == kInfinity)) {
      dof_output.inoperative_begin_execution_time = kInfinity;
      OverrideProfileTrace(dof_output.step1b);
    }

    if ((dof_output.inoperative_begin_execution_time != kInfinity) &&
        (dof_output.inoperative_end_execution_time != kInfinity)) {
      dof_output.inoperative_begin_execution_time -= kStep1TimeEpsilon;
      dof_output.inoperative_end_execution_time += kStep1TimeEpsilon;
    }

    dof_output.inoperative_begin_execution_time =
        std::max(dof_output.inoperative_begin_execution_time,
                 dof_output.min_execution_time);

    if (dof_output.inoperative_end_execution_time <=
        dof_output.inoperative_begin_execution_time) {
      double mean_inoperative_time =
          0.5 * (dof_output.inoperative_begin_execution_time +
                 dof_output.inoperative_end_execution_time);

      dof_output.inoperative_begin_execution_time =
          mean_inoperative_time - kStep1TimeEpsilon;
      dof_output.inoperative_end_execution_time =
          mean_inoperative_time + kStep1TimeEpsilon;

      if (dof_output.inoperative_begin_execution_time <
          dof_output.min_execution_time) {
        dof_output.inoperative_begin_execution_time = kInfinity;
        dof_output.inoperative_end_execution_time = kInfinity;

        OverrideProfileTrace(dof_output.step1b);
        OverrideProfileTrace(dof_output.step1c);
      }
    }

    if (dof_output.inoperative_end_execution_time <=
        dof_output.min_execution_time) {
      dof_output.inoperative_begin_execution_time = kInfinity;
      dof_output.inoperative_end_execution_time = kInfinity;

      OverrideProfileTrace(dof_output.step1b);
      OverrideProfileTrace(dof_output.step1c);
    }
  }
}

double CalculateExecutionTimeVector(const PositionInputs& inputs,
                                    Flags::SyncBehavior sync_behavior,
                                    Outputs& outputs) {
  double sync_time = GetMaximalMinExecutionTime(outputs);
  if (outputs.IsWithinAnInoperativeTimeInterval(sync_time)) {
    FixedVector<double, kMaxDofs> end_limits;
    for (const Outputs::DOF& dof : outputs.GetDOFs()) {
      if (dof.inoperative_end_execution_time > sync_time &&
          !outputs.IsWithinAnInoperativeTimeInterval(
              dof.inoperative_end_execution_time)) {
        end_limits.push_back(dof.inoperative_end_execution_time);
      }
    }
    sync_time = *absl::c_min_element(end_limits);
  }

  double user_defined_min_sync_time = inputs.GetMinimumSynchronizationTime();
  if (sync_time <= user_defined_min_sync_time) {
    // Set execution times for all dofs to the user time, as it must be at least
    // that.
    outputs.SetExecutionTimeIfSelected(user_defined_min_sync_time);
    return user_defined_min_sync_time;
  }

  // Set everything to the synchronization time if we're allowed to.
  if (sync_behavior != Flags::SyncBehavior::kNoSynchronization) {
    outputs.SetExecutionTimeIfSelected(sync_time);
    return sync_time;
  }

  // We're in async mode, so at least make sure we do the best we can with
  // respect to user defined min sync time.
  for (Outputs::DOF& dof : outputs.GetDOFs()) {
    double exec_time = dof.min_execution_time;
    if (exec_time < user_defined_min_sync_time) {
      if (dof.IsWithinAnInoperativeTimeInterval(user_defined_min_sync_time)) {
        exec_time = dof.inoperative_end_execution_time;
      } else {
        exec_time = user_defined_min_sync_time;
      }
    }
    dof.SetExecutionTimeIfSelected(exec_time);
  }

  return GetMaximalMinExecutionTime(outputs);
}
}  // namespace

bool CalcMinExecutionTime(const PositionInputs& inputs,
                          const PositionFlags& flags,
                          PositionOutputs& outputs) {
  // We track an array of whether or not step1a flipped
  MaxDOFFixedVector<bool> step1a_flipped(inputs.GetNumberOfDOFs());
  if (!CalcMinExecutionTimePerDOF(inputs, flags, outputs, step1a_flipped)) {
    return false;
  }

  // check if phase synchronization is expected
  // and if phase sync is possible
  if ((flags.synchronization_behavior ==
       Flags::SyncBehavior::kOnlyPhaseSynchronization) ||
      (flags.synchronization_behavior ==
       Flags::SyncBehavior::kPhaseSynchronizationIfPossible) ||
      (flags.synchronization_behavior ==
       Flags::SyncBehavior::kPhaseSynchronizationWhenCollinear)) {
    // run phase sync check
    // and if possible, compute phase sync data
    ComputePhaseSynchronizationData(true, inputs, outputs);
  } else {
    outputs.DisablePhaseSync();
  }

  // Calc inoperative zone irrespective of whether phase sync is possible or not
  if (!CalcInoperativeZoneStart(inputs, flags, step1a_flipped, outputs)) {
    return false;
  }

  if (!CalcInoperativeZoneEnd(inputs, flags, outputs)) {
    return false;
  }

  // Fixup the execution times computed above
  FixExecutionTimes(outputs);

  // Calculate synchronization time and execution time vector.
  outputs.SetSyncTime(CalculateExecutionTimeVector(
      inputs, flags.synchronization_behavior, outputs));

  return true;
}

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic
