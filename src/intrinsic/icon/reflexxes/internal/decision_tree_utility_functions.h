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

// This file contains common functions used to build decision trees.

#ifndef INTRINSIC_ICON_REFLEXXES_INTERNAL_DECISION_TREE_UTILITY_FUNCTIONS_H_
#define INTRINSIC_ICON_REFLEXXES_INTERNAL_DECISION_TREE_UTILITY_FUNCTIONS_H_

#include "intrinsic/icon/reflexxes/internal/motion_state.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

// Represents a decision in a decision tree.  If the decision is true, the
// true_func is called and passed in the given state; if it is false, false_func
// is called. In all cases, decision traces are applied using decision_num.
// Both true_func and _false_func must evaluate to callable types supporting
// StateT(const StateT&).
template <typename TrueFunc, typename FalseFunc, typename StateT>
StateT MakeDecision(const int decision_num, const bool decision,
                    const TrueFunc& true_func, const FalseFunc& false_func,
                    const StateT& state) {
  if (decision) {
    return AppendDecisionTrace(decision_num, true, state) | true_func;
  }
  return AppendDecisionTrace(decision_num, false, state) | false_func;
}

// Same as MakeDecision above, but instead of being passed a bool, a
// decision_func is passed and evaluated with the given state.  decision_func
// must be callable with the signature bool(const MotionState&).
template <typename DecisionFunc, typename TrueFunc, typename FalseFunc,
          typename StateT>
StateT MakeDecision(const int decision_num, const DecisionFunc& decision_func,
                    const TrueFunc& true_func, const FalseFunc& false_func,
                    const StateT& state) {
  return MakeDecision(decision_num, decision_func(state), true_func, false_func,
                      state);
}

// Makes a decision on whether or not applying the profile function profile_func
// results in a state that is in error.  If the state is not in error, that
// state is returned.  If the state is in error, error_func is applied to the
// original state, and that resulting state is returned.  Both profile_func and
// error_func must conform to callable objects with the signature StateT(const
// &StateT).
template <typename ProfileFunc, typename ErrorFunc, typename StateT>
StateT TryApplyProfileFunction(const int decision_num,
                               const ProfileFunc& profile_func,
                               const ErrorFunc& error_func,
                               const StateT& state) {
  auto s1 = profile_func(state);
  if (IsSuccess(s1)) {
    return AppendDecisionTrace(decision_num, true, s1);
  }
  return AppendDecisionTrace(decision_num, false, state) | error_func;
}

// Finishes a state successfully with infinite time.
constexpr auto FinishInfinite = SetT(kInfinity) | SetSuccess;

// Finishes a state in failure with infinite time.
constexpr auto FinishInfiniteFailure = SetT(kInfinity) | SetFailure;

// If the state is in failure, finish with infinite time.
constexpr auto IfErrorFinishInfinite = [](const auto& state) {
  if (IsFailure(state)) {
    return FinishInfinite(state);
  }
  return state;
};

// Writes the various states of the field back to the output substep struct.
template <typename State>
void WriteStateToSubStepOutput(const State& state,
                               Outputs::DOF::SubStep& output) {
  output.step_exec_time = GetT(state);
  output.applied_profile = GetProfile(state);
  output.result = true;
  output.decision_trace_size = GetDecisionTraceIndex(state);
  output.applied_profile_trace_size = GetDecisionProfileIndex(state);
}

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_INTERNAL_DECISION_TREE_UTILITY_FUNCTIONS_H_
