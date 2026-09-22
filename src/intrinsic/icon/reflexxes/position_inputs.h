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

#ifndef INTRINSIC_ICON_REFLEXXES_POSITION_INPUTS_H_
#define INTRINSIC_ICON_REFLEXXES_POSITION_INPUTS_H_

#include <utility>

#include "intrinsic/icon/reflexxes/inputs.h"

namespace intrinsic {
namespace reflexxes {

// Class for the input parameters of the position-based Online Trajectory
// Generation algorithm.
class PositionInputs : public Inputs {
 public:
  PositionInputs(const int num_dofs, const double cycle_time)
      : Inputs(num_dofs, cycle_time) {}

  // Useful for setting all values of the class.
  PositionInputs(const MaxDOFFixedVector<Inputs::DOF>& dofs,
                 const double cycle_time, const double min_sync_time)
      : Inputs(dofs, cycle_time, min_sync_time) {}

  std::pair<ErrorCodeForInvalidInputValues, int> CheckForValidity()
      const override;

  // Deselects all selected dofs that have invalid constraints.
  void DeselectInvalidDofs();

  // Checks whether the constraints are valid for all degrees of freedom.  Only
  // DOFs that are selected are considered for validity.
  bool CheckValidityOfConstraints() const;

  // Scale the input parameters, if necessary.
  // In order to ensure numerical stability, this method checks
  // whether the order of magnitude of the input parameters for each degree of
  // freedom is within  a predefined range of magnitude (c.f.
  // kLowerScalingThreshold and kUpperScalingThreshold), and scales them to the
  // predefined range, if they are not.
  // If one or more input values are out the range, in which numerical
  // stability is guaranteed, this function computes a scaling vector all
  // respective values.
  // Returns the scales used for each DOF.
  MaxDOFFixedVector<double> ScaleIfNecessary();
};

}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_POSITION_INPUTS_H_
