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

#ifndef INTRINSIC_ICON_REFLEXXES_INTERNAL_POSITION_H_
#define INTRINSIC_ICON_REFLEXXES_INTERNAL_POSITION_H_

#include "intrinsic/icon/reflexxes/internal/constants.h"
#include "intrinsic/icon/reflexxes/position_flags.h"
#include "intrinsic/icon/reflexxes/position_inputs.h"
#include "intrinsic/icon/reflexxes/position_outputs.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

struct PositionState {
  PositionState(const int num_dofs, const double cycle_time)
      : inputs(num_dofs, cycle_time), outputs(num_dofs, cycle_time) {}

  PositionInputs inputs;
  PositionFlags flags;
  PositionOutputs outputs;
};

// Calculates the minimum execution time for each DOF using the Position OTG
// algorithm.
bool CalcMinExecutionTimePerDOF(const PositionInputs& inputs,
                                const PositionFlags& flags,
                                PositionOutputs& outputs,
                                MaxDOFFixedVector<bool>& step1a_flipped);

// Calculates the start time of the inoperative zone (if any) for each DOF using
// the Position OTG algorithm.
bool CalcInoperativeZoneStart(const PositionInputs& inputs,
                              const PositionFlags& flags,
                              const MaxDOFFixedVector<bool>& step1a_flipped,
                              PositionOutputs& outputs);

// Calculates the end time of the inoperative zone (if any) for each DOF using
// the Position OTG algorithm.
bool CalcInoperativeZoneEnd(const PositionInputs& inputs,
                            const PositionFlags& flags,
                            PositionOutputs& outputs);

// Calculates the minimum execution time across for the Position OTG algorithm
// for all DOFs.
bool CalcMinExecutionTime(const PositionInputs& inputs,
                          const PositionFlags& flags, PositionOutputs& outputs);

// Synchronizes the DOFs and calculates the motion trajectories for the Position
// OTG algorithm.
bool SynchronizeDOFs(const Inputs& inputs, const PositionFlags& flags,
                     PositionOutputs& outputs);

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_INTERNAL_POSITION_H_
