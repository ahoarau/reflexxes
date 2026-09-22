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

#include "intrinsic/icon/reflexxes/internal/calc_min_execution_time_utils.h"

#include <cmath>
#include <cstdlib>

#include "absl/algorithm/container.h"
#include "intrinsic/icon/reflexxes/inputs.h"
#include "intrinsic/icon/reflexxes/internal/constants.h"
#include "intrinsic/icon/reflexxes/internal/math.h"
#include "intrinsic/icon/reflexxes/outputs.h"
#include "intrinsic/util/fixed_vector.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

namespace {

// Get the number of selected degrees of freedom
int GetNumberOfSelectedDOFs(const Inputs& inputs) {
  return absl::c_count_if(inputs.GetDOFs(),
                          [](const Inputs::DOF& dof) { return dof.selected; });
}

auto GetMaximalMinExecutionTimeDOF(const Outputs& outputs) {
  auto dofs = outputs.GetDOFs();
  return absl::c_max_element(
      dofs, [](const Outputs::DOF& a, const Outputs::DOF& b) {
        return a.min_execution_time < b.min_execution_time;
      });
}

int GetMaximalMinExecutionTimeDofIndex(const Outputs& outputs) {
  if (outputs.GetNumberOfDOFs() == 0) {
    return 0;
  }
  return GetMaximalMinExecutionTimeDOF(outputs)->index;
}

bool Normalize(MaxDOFFixedVector<double>& values) {
  double magnitude = 0.0;
  for (double value : values) {
    magnitude += Power2(value);
  }
  magnitude = GetSqrt(magnitude);
  bool has_magnitude = magnitude >= kAbsPhaseSyncEpsilon;

  if (has_magnitude) {
    for (double& value : values) {
      value /= magnitude;
    }
  }

  return has_magnitude;
}

}  // namespace

double GetMaximalMinExecutionTime(const Outputs& outputs) {
  if (outputs.GetNumberOfDOFs() == 0) {
    return 0.0;
  }
  return GetMaximalMinExecutionTimeDOF(outputs)->min_execution_time;
}

void ComputePhaseSynchronizationData(const bool consider_delta_pos,
                                     const Inputs& inputs, Outputs& outputs) {
  // Check if there is only one selected dof. If yes phase sync is trivially
  // possible.
  if (GetNumberOfSelectedDOFs(inputs) == 1) {
    auto dofs = outputs.GetDOFs();
    auto dof_iter = absl::c_find_if(
        dofs, [](const Outputs::DOF& dof) { return dof.selected; });
    dof_iter->phase_sync_scale = 1.0;
    outputs.EnablePhaseSync(dof_iter->index);
    return;
  }

  // We track four vectors, each over all DOFs for a specific attribute: delta
  // position, velocity, acceleration, and target velocity.
  FixedVector<MaxDOFFixedVector<double>, 4> normalized_values;

  MaxDOFFixedVector<double> values(inputs.GetNumberOfDOFs());
  for (const Inputs::DOF& dof : inputs.GetDOFs()) {
    values[dof.index] = dof.selected && consider_delta_pos
                            ? dof.target_position - dof.position
                            : 0.;
  }
  if (Normalize(values)) {
    normalized_values.emplace_back(values);
  }

  for (const Inputs::DOF& dof : inputs.GetDOFs()) {
    values[dof.index] = dof.selected ? dof.velocity : 0.;
  }
  if (Normalize(values)) {
    normalized_values.emplace_back(values);
  }

  for (const Inputs::DOF& dof : inputs.GetDOFs()) {
    values[dof.index] = dof.selected ? dof.acceleration : 0.;
  }
  if (Normalize(values)) {
    normalized_values.emplace_back(values);
  }

  for (const Inputs::DOF& dof : inputs.GetDOFs()) {
    values[dof.index] = dof.selected ? dof.target_velocity : 0.;
  }
  if (Normalize(values)) {
    normalized_values.emplace_back(values);
  }

  // If there are no vectors with magnitude, exit early.
  if (normalized_values.empty()) {
    if (consider_delta_pos) {
      // Set the default phase sync dof and scale of unity for all dofs.
      for (Outputs::DOF& output_dof : outputs.GetDOFs()) {
        output_dof.phase_sync_scale = 1.0;
      }

      outputs.EnablePhaseSync(0);
      return;
    }

    // In velocity mode this does not let us phase sync.
    outputs.DisablePhaseSync();
    return;
  }

  // Check that all vectors are collinear with the first one.  (Its
  // transitive, so we don't need to test everything).
  for (auto iter = normalized_values.begin() + 1;
       iter != normalized_values.end(); ++iter) {
    // If the inner product (dot) is 1.0, they are collinear.
    if (!EpsilonEqual(
            fabs(absl::c_inner_product(normalized_values.front(), *iter, 0.0)),
            1.0, kAbsPhaseSyncEpsilon)) {
      // If anything is not collinear, disable and bail.
      outputs.DisablePhaseSync();
      return;
    }
  }

  // Find the phase sync dof, which is the dof with the largest min execution
  // time.
  int phase_sync_dof = GetMaximalMinExecutionTimeDofIndex(outputs);
  outputs.EnablePhaseSync(phase_sync_dof);

  // Find the reference value used for scale factor computation.
  auto reference_iter = absl::c_find_if(
      normalized_values, [&phase_sync_dof](const auto& normalized_value) {
        return fabs(normalized_value[phase_sync_dof]) >= kAbsPhaseSyncEpsilon;
      });

  // Calculate the scale factor for each dof.
  if (reference_iter != normalized_values.end()) {
    for (Outputs::DOF& dof : outputs.GetDOFs()) {
      if (dof.selected && (dof.index != phase_sync_dof)) {
        dof.phase_sync_scale =
            (*reference_iter)[dof.index] / (*reference_iter)[phase_sync_dof];
      } else {
        dof.phase_sync_scale = 1.0;
      }
    }
  }
}

void ScaleLimitsForPhaseSync(Inputs& inputs, const Outputs& outputs) {
  // The index of the DoF for which the highest phase sync scaling
  // factor has been computed.
  int max_scale_factor_dof_idx = -1;
  // The maximum scale factor from all DoFs, to be computed below - if it
  // exists.
  double max_scale_factor = 0.0;

  for (const Inputs::DOF& dof : inputs.GetDOFs()) {
    // If applicable, update the max scale factor and the current lead DoF
    // index.
    if (dof.selected &&
        (std::abs(outputs.GetDOFs()[dof.index].phase_sync_scale) >
         max_scale_factor)) {
      max_scale_factor =
          std::abs(outputs.GetDOFs()[dof.index].phase_sync_scale);
      max_scale_factor_dof_idx = dof.index;
    }
  }

  // Return early in case there is no admissible scaling to be applied.
  if (max_scale_factor_dof_idx < 0) return;

  // For every DoF affected by phase synchronization, scale limits according
  // to the above determined leading DoF, compare Eq. 6.14 in "On-Line
  // Trajectory Generation in Robotic Systems" by T. Kroeger. Note that only
  // "downscaling" limits is appropriate, this is why we do not employ the
  // 'phase_sync_dof' from 'outputs', as other code portions do.
  const auto& lead_dof = inputs.GetDOFs()[max_scale_factor_dof_idx];
  for (Inputs::DOF& dof : inputs.GetDOFs()) {
    if (dof.selected && dof.index != max_scale_factor_dof_idx) {
      const double limit_scale_factor = std::abs(
          outputs.GetDOFs()[dof.index].phase_sync_scale / max_scale_factor);

      dof.min_velocity = limit_scale_factor * lead_dof.min_velocity;
      dof.max_velocity = limit_scale_factor * lead_dof.max_velocity;
      dof.min_acceleration = limit_scale_factor * lead_dof.min_acceleration;
      dof.max_acceleration = limit_scale_factor * lead_dof.max_acceleration;
      dof.min_jerk = limit_scale_factor * lead_dof.min_jerk;
      dof.max_jerk = limit_scale_factor * lead_dof.max_jerk;
    }
  }
}

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic
