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

#include "intrinsic/icon/reflexxes/inputs.h"

#include "absl/log/check.h"

namespace intrinsic {
namespace reflexxes {

void Inputs::DOF::ScaleValues(const double scale_factor) {
  CHECK(scale_factor > 0.0);

  position = position * scale_factor;
  velocity = velocity * scale_factor;
  acceleration = acceleration * scale_factor;
  // bounds
  max_position = max_position * scale_factor;
  max_velocity = max_velocity * scale_factor;
  max_acceleration = max_acceleration * scale_factor;
  max_jerk = max_jerk * scale_factor;
  min_position = min_position * scale_factor;
  min_velocity = min_velocity * scale_factor;
  min_acceleration = min_acceleration * scale_factor;
  min_jerk = min_jerk * scale_factor;
  // target state
  target_position = target_position * scale_factor;
  target_velocity = target_velocity * scale_factor;
  alt_target_velocity = alt_target_velocity * scale_factor;
}

Inputs::Inputs(const int num_dofs, const double cycle_time)
    : dofs_(num_dofs), cycle_time_(cycle_time), min_sync_time_(0.0) {
  for (int i = 0; i < dofs_.size(); ++i) {
    dofs_[i].index = i;
  }
}

Inputs::Inputs(const MaxDOFFixedVector<DOF>& dofs, const double cycle_time,
               const double min_sync_time)
    : dofs_(dofs), cycle_time_(cycle_time), min_sync_time_(min_sync_time) {
  for (int i = 0; i < dofs_.size(); ++i) {
    dofs_[i].index = i;
  }
}

Inputs::DOF FlipInputParameters(const Inputs::DOF& i) {
  Inputs::DOF result = i;
  result.position *= -1;
  result.velocity *= -1;
  result.acceleration *= -1;
  result.target_position *= -1;
  result.target_velocity *= -1;
  result.max_velocity = -i.min_velocity;
  result.min_velocity = -i.max_velocity;
  result.min_acceleration = -i.max_acceleration;
  result.max_acceleration = -i.min_acceleration;
  result.min_jerk = -i.max_jerk;
  result.max_jerk = -i.min_jerk;
  return result;
}

}  // namespace reflexxes
}  // namespace intrinsic
