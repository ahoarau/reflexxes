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

#ifndef INTRINSIC_ICON_REFLEXXES_VELOCITY_INPUTS_H_
#define INTRINSIC_ICON_REFLEXXES_VELOCITY_INPUTS_H_

#include <utility>

#include "intrinsic/icon/reflexxes/inputs.h"

namespace intrinsic {
namespace reflexxes {

// Class for the input parameters of the velocity-based Online Trajectory
// Generation algorithm
class VelocityInputs : public Inputs {
 public:
  VelocityInputs(const int num_dofs, const double cycle_time)
      : Inputs(num_dofs, cycle_time) {}

  explicit VelocityInputs(const Inputs& input_param) : Inputs(input_param) {}

  // Constructor to initialize all info.
  VelocityInputs(const MaxDOFFixedVector<Inputs::DOF>& dofs,
                 const double cycle_time, const double min_sync_time)
      : Inputs(dofs, cycle_time, min_sync_time) {}

  // Copy operator to copy over the base class inputs.
  VelocityInputs& operator=(const Inputs& input_param) {
    Inputs::operator=(input_param);
    return *this;
  }

  std::pair<ErrorCodeForInvalidInputValues, int> CheckForValidity()
      const override;

  static bool CheckValidityOfConstraintsForSingleDOF(const Inputs::DOF& dof);

  // Checks whether the constraints are valid for all selected DOFs.
  bool CheckValidityOfConstraints() const;

  // Sets the selection vector to true for all selected dofs that have valid
  // constraints, and false for those who do not.  If a dof is not already
  // selected there is no change.
  void DeselectInvalidDofs();

  // Scale the input parameters, if necessary.
  // In order to ensure numerical stability, this method checks
  // whether the order of magnitude of  the input parameters for each degree of
  // freedom is within a predefined range of magnitude (c.f.
  // kLowerScalingThreshold and kUpperScalingThreshold), and scales them to the
  // predefined range, if they are not. If one or more input values are out the
  // range, in which numerical stability is guaranteed, this function computes a
  // scaling vector all respective values. This vector is subsequently used in
  // TypeIVRMLVelocity, such that the  Trajectory Generation algorithm can
  // internally work with the scaled and linearly transformed values.
  MaxDOFFixedVector<double> ScaleIfNecessary();
};

}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_VELOCITY_INPUTS_H_
