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

#ifndef INTRINSIC_ICON_REFLEXXES_VELOCITY_OUTPUTS_H_
#define INTRINSIC_ICON_REFLEXXES_VELOCITY_OUTPUTS_H_

#include <cstdint>

#include "intrinsic/icon/reflexxes/outputs.h"
#include "intrinsic/icon/reflexxes/velocity_inputs.h"

namespace intrinsic {
namespace reflexxes {

// Class for the output parameters of the velocity-based Online Trajectory
// Generation algorithm
class VelocityOutputs : public Outputs {
 public:
  VelocityOutputs(const unsigned int num_dofs, const double cycle_time)
      : Outputs(num_dofs, cycle_time),
        position_values_at_target_velocity_(num_dofs) {}

  // Explicit "init everything" constructor.  Should not generally be used.
  VelocityOutputs(
      const MaxDOFFixedVector<DOF>& dofs, const double cycle_time,
      const double sync_time, const uint64_t time_counter, const Status status,
      const bool is_valid_output_available, const bool phase_sync_enabled,
      const int phase_sync_dof_index,
      const MaxDOFFixedVector<double>& position_values_at_target_velocity)
      : Outputs(dofs, cycle_time, sync_time, time_counter, status,
                is_valid_output_available, phase_sync_enabled,
                phase_sync_dof_index),
        position_values_at_target_velocity_(
            position_values_at_target_velocity) {}

  // Assign from the base class.
  VelocityOutputs& operator=(const Outputs& output_param) {
    Outputs::operator=(output_param);
    return *this;
  }

  // Gets the position value for each DOF at the instant when its target
  // velocity is reached
  MaxDOFFixedVector<double>& GetPositionValuesAtTargetVelocity() {
    return position_values_at_target_velocity_;
  }

  // Gets the position value for each DOF at the instant when its target
  // velocity is reached.
  const MaxDOFFixedVector<double>& GetPositionValuesAtTargetVelocity() const {
    return position_values_at_target_velocity_;
  }

 private:
  // the position value for each DOF at the instant when its target velocity is
  // reached
  MaxDOFFixedVector<double> position_values_at_target_velocity_;
};

}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_VELOCITY_OUTPUTS_H_
