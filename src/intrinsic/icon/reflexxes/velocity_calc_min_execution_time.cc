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

#include "intrinsic/icon/reflexxes/flags.h"
#include "intrinsic/icon/reflexxes/internal/calc_min_execution_time_utils.h"
#include "intrinsic/icon/reflexxes/internal/decision_tree_utility_functions.h"
#include "intrinsic/icon/reflexxes/internal/functional.h"
#include "intrinsic/icon/reflexxes/internal/intermediate_profile_funcs.h"
#include "intrinsic/icon/reflexxes/internal/motion_state.h"
#include "intrinsic/icon/reflexxes/internal/motion_state_transforms.h"
#include "intrinsic/icon/reflexxes/outputs.h"
#include "intrinsic/icon/reflexxes/profile.h"
#include "intrinsic/icon/reflexxes/velocity_flags.h"
#include "intrinsic/icon/reflexxes/velocity_inputs.h"
#include "intrinsic/icon/reflexxes/velocity_outputs.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

using MotionState = MotionStateT<>;

namespace {

// v->vtrgt (PosTrap)
constexpr auto CalcExecTimePosTrap =
    VUpToVTrgtPosTrap | SetProfile(Profile::kPosTrap);

// v->vtrgt (PosTri)
constexpr auto CalcExecTimePosTri =
    VUpToVTrgtPosTri | SetProfile(Profile::kPosTri);

// The decision tree.
// This is done as a struct to allow the decisions to be listed in order
// starting from the first decision.  Bundling the decisions into a struct
// allows us to refer to later decisions before they are declared, whereas a
// namespace wouldn't allow that.
struct Decisions {
  static inline MotionState D1(const MotionState& s) {
    return MakeDecision(1, CmpAGteZero, &D2, FlipState | &D2, s);
  }

  static inline MotionState D2(const MotionState& s) {
    return MakeDecision(2, CmpALteAMax, &D3,
                        CalcExecTimeNegLinAdowntoAMax | &D3, s);
  }

  static inline MotionState D3(const MotionState& s) {
    return MakeDecision(3, ADownToZero | CmpVLteVTrgt, &D4,
                        CalcExecTimeNegLinAdowntoZero | FlipState | &D4, s);
  }

  static inline MotionState D4(const MotionState& s) {
    return MakeDecision(4, AUpToAMaxDownToZero | CmpVLteVTrgt,
                        CalcExecTimePosTrap, CalcExecTimePosTri, s);
  }
};

}  // namespace

bool CalcMinExecutionTime(const VelocityInputs& inputs,
                          const VelocityFlags& flags,
                          VelocityOutputs& outputs) {
  for (const auto [dof_input, dof_output] :
       Zip(inputs.GetDOFs(), outputs.GetDOFs())) {
    // Check if the dof is available for use and if the constraints are valid.
    if (!dof_input.selected |
        !VelocityInputs::CheckValidityOfConstraintsForSingleDOF(dof_input)) {
      dof_output.step1a.result = true;

      continue;
    }

    MotionState s = Decisions::D1(InitializeMotionState(
        StateBase(dof_input, dof_output, dof_output.step1a, nullptr),
        dof_input));

    if (IsFailure(s)) {
      // Any error on any DOF is enough to exit the decision, no need to
      // continue.
      return false;
    }

    WriteStateToSubStepOutput(s, dof_output.step1a);
    dof_output.min_execution_time = GetT(s);
  }

  // Check if phase synchronization is expected and if phase sync is possible.
  if ((flags.synchronization_behavior ==
       Flags::SyncBehavior::kOnlyPhaseSynchronization) ||
      (flags.synchronization_behavior ==
       Flags::SyncBehavior::kPhaseSynchronizationIfPossible) ||
      (flags.synchronization_behavior ==
       Flags::SyncBehavior::kPhaseSynchronizationWhenCollinear)) {
    ComputePhaseSynchronizationData(false, inputs, outputs);
  } else {
    outputs.DisablePhaseSync();
  }

  const double user_defined_sync_time = inputs.GetMinimumSynchronizationTime();
  if (flags.synchronization_behavior ==
      Flags::SyncBehavior::kNoSynchronization) {
    // If there's no synchronization, ensure the execution times are at least
    // the requested minimum time.
    double synchronization_time = 0.0;
    for (Outputs::DOF& d : outputs.GetDOFs()) {
      d.execution_time = std::max(d.min_execution_time, user_defined_sync_time);
      synchronization_time = std::max(synchronization_time, d.execution_time);
    }

    outputs.SetSyncTime(synchronization_time);
    return true;
  }

  // Fill out the synchronized execution time and set it across all DOFs.
  double synchronization_time = std::max(
      outputs.IsPhaseSyncEnabled()
          ? outputs.GetDOFs()[outputs.GetPhaseSyncDOFIndex()].min_execution_time
          : GetMaximalMinExecutionTime(outputs),
      user_defined_sync_time);
  outputs.SetExecutionTimeIfSelected(synchronization_time);
  outputs.SetSyncTime(synchronization_time);

  return true;
}

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic
