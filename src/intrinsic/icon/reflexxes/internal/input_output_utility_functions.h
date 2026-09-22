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

// This file provides a collection of utilities that work with the Inputs and
// Outputs classes and their derivatives.

#ifndef INTRINSIC_ICON_REFLEXXES_INTERNAL_INPUT_OUTPUT_UTILITY_FUNCTIONS_H_
#define INTRINSIC_ICON_REFLEXXES_INTERNAL_INPUT_OUTPUT_UTILITY_FUNCTIONS_H_

#include "intrinsic/icon/reflexxes/inputs.h"
#include "intrinsic/icon/reflexxes/internal/constants.h"
#include "intrinsic/icon/reflexxes/outputs.h"
#include "intrinsic/icon/reflexxes/position_inputs.h"
#include "intrinsic/icon/reflexxes/position_outputs.h"
#include "intrinsic/icon/reflexxes/velocity_inputs.h"
#include "intrinsic/icon/reflexxes/velocity_outputs.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

// Checks if the PVA parameters from the inputs match those of the outputs.
bool DoOutputsDeviateFromInputs(const Inputs& inputs, const Outputs& outputs);

// Returns true if the limits, targets, and timing information inthe old_inputs
// are compatible with the new_inputs.
bool AreInputsCompatible(const Inputs& old_inputs, const Inputs& new_inputs);

// Returns true if the limits, targets, and timing information inthe old_inputs
// are compatible with the new_inputs.
bool AreInputsCompatible(const PositionInputs& old_inputs,
                         const PositionInputs& new_inputs);

// Checks whether the output trajectory breaches the positional limits. This
// method is not called in the algorithmic core of Reflexxes, it is only used
// for checking if a fallback-behaviour stop needs to be requested. To allow for
// robust checks despite small inaccuracies during problem scaling, positional
// checks need to have a tolerance.
bool ArePositionalLimitsBreached(
    const Inputs& inputs, const Outputs& outputs,
    double tolerance = kPositionLimitBreachedThreshold);

// Checks whether the output trajectory breaches the velocity limits
bool AreVelocityLimitsBreached(const Inputs& inputs, const Outputs& outputs);

// Checks whether the output trajectory breaches the acceleration limits
bool AreAccelerationLimitsBreached(const Inputs& inputs,
                                   const Outputs& outputs);

// Checks whether the output trajectory breaches the jerk limits
bool AreJerkLimitsBreached(const Inputs& inputs, const Outputs& outputs);

// Computes minimum and maximum positions for each dof, their associated
// times and the states of other degrees of freedom at these times.
void ComputeExtremaData(const PositionInputs& inputs, PositionOutputs& outputs);
void ComputeExtremaData(const VelocityInputs&, VelocityOutputs& outputs);

// Copies the selection vector from inputs to outputs.
void CopySelectionVector(const Inputs& inputs, Outputs& outputs);

// Copies the new state of the output to the input.
void CopyNewStateToInputs(const Outputs& outputs, Inputs& inputs);

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_INTERNAL_INPUT_OUTPUT_UTILITY_FUNCTIONS_H_
