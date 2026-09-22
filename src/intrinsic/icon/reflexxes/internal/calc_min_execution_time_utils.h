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

#ifndef INTRINSIC_ICON_REFLEXXES_INTERNAL_CALC_MIN_EXECUTION_TIME_UTILS_H_
#define INTRINSIC_ICON_REFLEXXES_INTERNAL_CALC_MIN_EXECUTION_TIME_UTILS_H_

#include "intrinsic/icon/reflexxes/inputs.h"
#include "intrinsic/icon/reflexxes/outputs.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

// Returns the largest min execution time over all the dofs.
double GetMaximalMinExecutionTime(const Outputs& outputs);

// Computes whether or not phase sync is enabled, and stores that in the output
// parameters.
void ComputePhaseSynchronizationData(bool consider_delta_pos,
                                     const Inputs& inputs, Outputs& outputs);

// Scales individual DoF limits such that they match the scaling factors for
// phase synchronization. Scaled limits are a function of the limits of the DoF
// which "leads" the phase synchronization. This is required to achieve
// consistent trajectories, see Chapter 6 in "On-Line Trajectory Generation in
// Robotic Systems" by T. Kroeger. Returns 'true' if limit scaling was
// successful, 'false' otherwise.
void ScaleLimitsForPhaseSync(Inputs& inputs, const Outputs& outputs);

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_INTERNAL_CALC_MIN_EXECUTION_TIME_UTILS_H_
