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

#ifndef INTRINSIC_ICON_REFLEXXES_POSITION_OUTPUTS_H_
#define INTRINSIC_ICON_REFLEXXES_POSITION_OUTPUTS_H_

#include <cstdint>

#include "intrinsic/icon/reflexxes/outputs.h"
#include "intrinsic/icon/reflexxes/position_flags.h"
#include "intrinsic/icon/reflexxes/position_inputs.h"

namespace intrinsic {
namespace reflexxes {

// Class for the output parameters of the position-based Online Trajectory
// Generation algorithm.
class PositionOutputs : public Outputs {
 public:
  PositionOutputs(const unsigned int num_dofs, const double cycle_time)
      : Outputs(num_dofs, cycle_time) {}

  PositionOutputs(const MaxDOFFixedVector<DOF>& dofs, const double cycle_time,
                  const double sync_time, const uint64_t time_counter,
                  const Status status, const bool is_valid_output_available,
                  const bool phase_sync_enabled, const int phase_sync_dof_index,
                  const bool trajectory_exceeds_target_position)
      : Outputs(dofs, cycle_time, sync_time, time_counter, status,
                is_valid_output_available, phase_sync_enabled,
                phase_sync_dof_index),
        trajectory_exceeds_target_position_(
            trajectory_exceeds_target_position) {}

  // Copy operator to pull data in from the base class.
  PositionOutputs& operator=(const Outputs& output_param) {
    Outputs::operator=(output_param);
    SetTrajectoryExceedsTargetPosition(false);
    return *this;
  }

  // Returns the boolean flags that indicates whether the currently computed
  // trajectory will exceed the desired target position in order to reach the
  // desired target state.
  bool GetTrajectoryExceedsTargetPosition() const {
    return trajectory_exceeds_target_position_;
  }

  // Sets the boolean flags that indicates, whether the currently computed
  // trajectory will exceed the desired target position in order to reach the
  // desired target state.
  void SetTrajectoryExceedsTargetPosition(const bool val) {
    trajectory_exceeds_target_position_ = val;
  }

 private:
  // Indicates, whether the currently computed trajectory will exceed the
  // desired target position in order to reach the desired target state.
  bool trajectory_exceeds_target_position_ = false;
};

}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_POSITION_OUTPUTS_H_
