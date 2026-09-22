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

#ifndef INTRINSIC_ICON_REFLEXXES_INTERNAL_POSITION_CALC_MIN_EXECUTION_TIME_BASE_H_
#define INTRINSIC_ICON_REFLEXXES_INTERNAL_POSITION_CALC_MIN_EXECUTION_TIME_BASE_H_

#include <functional>

#include "absl/strings/str_format.h"
#include "intrinsic/icon/reflexxes/inputs.h"
#include "intrinsic/icon/reflexxes/internal/motion_state.h"
#include "intrinsic/icon/reflexxes/internal/motion_state_transforms.h"
#include "intrinsic/icon/reflexxes/outputs.h"
#include "intrinsic/icon/reflexxes/position_flags.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

using MotionState = MotionStateT<const bool>;
enum PositionStep1MotionStateMembers {
  kStep1AFlipped = details::MotionStateMembers::kMaxBaseProperties
};

constexpr auto GetStep1AProfile = [](const MotionState& s) {
  return GetDOFOutput(s).step1a.applied_profile;
};

constexpr auto GetStep1AAppliedProfileFlipped = [](const MotionState& s) {
  return std::get<PositionStep1MotionStateMembers::kStep1AFlipped>(s);
};

constexpr auto SetStep1AAppliedProfileFlipped =
    Curry<bool, MotionState>([](const bool flipped, const MotionState& s) {
      return details::SetValue<PositionStep1MotionStateMembers::kStep1AFlipped>(
          flipped, s);
    });

// Helper to create the state for a position step 1 decision tree and run the
// given tree via the passed in decsion_tree_func.
inline MotionState RunPositionStep1DecisionTree(
    const std::function<MotionState(const MotionState&)>& decision_tree_func,
    const Inputs::DOF& dof_input, const PositionFlags& flags,
    Outputs::DOF& dof_output, Outputs::DOF::SubStep& sub_step,
    const bool step1a_flipped = false) {
  StateBase state_base(
      dof_input, dof_output, sub_step, nullptr, 0,
      flags.behavior_if_initial_state_breaches_constraints ==
          PositionFlags::BehaviorIfInitialStateBreachesConstraints::
              kGetIntoBoundariesFast);
  return decision_tree_func(
      InitializeMotionState<const bool>(state_base, dof_input, step1a_flipped));
}

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_INTERNAL_POSITION_CALC_MIN_EXECUTION_TIME_BASE_H_
