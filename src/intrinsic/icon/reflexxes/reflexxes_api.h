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

// The primary interface to the Reflexxes Online Trajectory Generation
// algorithm.
//
// Typically, usage follows this pattern:
//
// 1. Setup input parameters.  This includes filling out the current state in
// the input parameters (see inputs.h, position_inputs.h, velocity_inputs.h) and
// the desired option flags (see flags.h, position_flags.h, veclocity_flags.h)
// for the call.  All parameter objects, including inputs, outputs, and the
// state object are safe to create on the stack or can be heap allocated prior
// to a realtime loop.  The State object must also be created and passed in.
//
// 2. Call the ComputePosition or ComputeVelocity function and evaluate the
// returned Status for errors (see status.h).  A success call is either kWorking
// (still something to do) or kFinalStateReached (done).  Failure statuses and
// recovery depend on the types of flags set in the passed in Flags structure.
//
// 3. If the call was successful, the outputs (see outputs.h,
// position_outputs.h, velocity_outputs.h) can be queried for the current
// estimated state for each degree of freedom.  Queries can also be made about
// the overall plan to reach target, including motion profiles for each segment
// of the plan.
//
// 4. The loop begins again, and the previous (or actual) state is fed back into
// the input.

#ifndef INTRINSIC_ICON_REFLEXXES_REFLEXXES_API_H_
#define INTRINSIC_ICON_REFLEXXES_REFLEXXES_API_H_

#include "intrinsic/icon/reflexxes/internal/constants.h"
#include "intrinsic/icon/reflexxes/internal/position.h"
#include "intrinsic/icon/reflexxes/internal/velocity.h"
#include "intrinsic/icon/reflexxes/position_flags.h"
#include "intrinsic/icon/reflexxes/position_inputs.h"
#include "intrinsic/icon/reflexxes/position_outputs.h"
#include "intrinsic/icon/reflexxes/status.h"
#include "intrinsic/icon/reflexxes/velocity_flags.h"
#include "intrinsic/icon/reflexxes/velocity_inputs.h"
#include "intrinsic/icon/reflexxes/velocity_outputs.h"

namespace intrinsic {
namespace reflexxes {

// A state object to be passed to one of the Compute calls below.
// The contents of this structure should be considered opaque and not relied
// upon in any way.
struct State {
  State(int num_dofs, double cycle_time);

  int GetNumberOfDOFs() const;
  double GetCycleTime() const;

 private:
  friend Status ComputePosition(const PositionInputs&, const PositionFlags&,
                                PositionOutputs&, State&);
  friend Status ComputeVelocity(const VelocityInputs&, const VelocityFlags&,
                                VelocityOutputs&, State&);

  internal::PositionState position_state;
  internal::VelocityState velocity_state;
  internal::VelocityState fallback_velocity_state;
};

// The position-based Type IV Online Trajectory Generation algorithm.
// Computes the result of applying the algorithm with the given inputs and flags
// and stores the result in the output parameter.
// The state param is an in/out variable that is passed to the initial call of
// ComputePosition, computed, and then passed to subsequent calls to update the
// same trajectory.  This is an opaque structure. After initial construction,
// the library will take care of any re-initialization internally.
Status ComputePosition(const PositionInputs& inputs, const PositionFlags& flags,
                       PositionOutputs& outputs, State& state);

// The velocity-based Type IV Online Trajectory Generation algorithm.
// Computes the result of applying the algorithm with the given inputs and flags
// and stores the result in the output parameter.
// The state param is an in/out variable that is passed to the initial call of
// ComputeVelocity, computed, and then passed to subsequent calls to update the
// same trajectory.  This is an opaque structure. After initial construction,
// the library will take care of any re-initialization internally.
Status ComputeVelocity(const VelocityInputs& inputs, const VelocityFlags& flags,
                       VelocityOutputs& outputs, State& state);

// Exposes whether a stop trajectory for the given inputs terminates in a
// state that respects positional limits.
bool DoesStopMotionFromNewStateBreachPositionalLimits(
    const VelocityInputs& inputs, const VelocityFlags& flags,
    const VelocityOutputs& outputs);

}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_REFLEXXES_API_H_
