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

#ifndef INTRINSIC_ICON_REFLEXXES_INTERNAL_VELOCITY_H_
#define INTRINSIC_ICON_REFLEXXES_INTERNAL_VELOCITY_H_

#include "intrinsic/icon/reflexxes/internal/constants.h"
#include "intrinsic/icon/reflexxes/velocity_flags.h"
#include "intrinsic/icon/reflexxes/velocity_inputs.h"
#include "intrinsic/icon/reflexxes/velocity_outputs.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

struct VelocityState {
  VelocityState(const int num_dofs, const double cycle_time)
      : inputs(num_dofs, cycle_time), outputs(num_dofs, cycle_time) {}

  VelocityInputs inputs;
  VelocityFlags flags;
  VelocityOutputs outputs;
};

// Calculates the minimum execution time across for the Position OTG algorithm
// for all DOFs.
bool CalcMinExecutionTime(const VelocityInputs& inputs,
                          const VelocityFlags& flags, VelocityOutputs& outputs);

// Synchronizes the DOFs and calculates the motion trajectories for the Position
// OTG algorithm.
bool SynchronizeDOFs(const Inputs& inputs, const VelocityFlags& flags,
                     VelocityOutputs& outputs);

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_INTERNAL_VELOCITY_H_
