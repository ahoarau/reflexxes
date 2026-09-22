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

// This file defines the decision tree for the CalcInoperativeZoneStart
// function, otherwise known as "Position Step1B" in the Reflexxes literature.
// See
// https://github.com/intrinsic-ai/intrinsic-core/blob/main/intrinsic_control/intrinsic/icon/reflexxes/g3doc/step1b_decision_tree.pdf
// for a detailed description of what's happening in this file.

#include <cmath>
#include <utility>

#include "intrinsic/icon/reflexxes/inputs.h"
#include "intrinsic/icon/reflexxes/internal/constants.h"
#include "intrinsic/icon/reflexxes/internal/decision_tree_utility_functions.h"
#include "intrinsic/icon/reflexxes/internal/functional.h"
#include "intrinsic/icon/reflexxes/internal/intermediate_profile_funcs.h"
#include "intrinsic/icon/reflexxes/internal/math.h"
#include "intrinsic/icon/reflexxes/internal/motion_state.h"
#include "intrinsic/icon/reflexxes/internal/motion_state_transforms.h"
#include "intrinsic/icon/reflexxes/internal/position.h"
#include "intrinsic/icon/reflexxes/internal/position_calc_exec_time_funcs.h"
#include "intrinsic/icon/reflexxes/internal/position_calc_min_execution_time_base.h"
#include "intrinsic/icon/reflexxes/position_flags.h"
#include "intrinsic/icon/reflexxes/position_inputs.h"
#include "intrinsic/icon/reflexxes/position_outputs.h"
#include "intrinsic/icon/reflexxes/profile.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

namespace {

Range UpdateLimitsForPosTriFunc(const double execution_time,
                                const Range& param_min_max,
                                const MotionState& s) {
  auto p_func = GetPositionErrorForPosTriNegTri(s);
  auto dp_func = GetDerivativePositionErrorForPosTriNegTri(s);
  double derivative_position_error_at_param_min =
      dp_func(param_min_max.min, nullptr, nullptr);
  double derivative_position_error_at_param_max =
      dp_func(param_min_max.max, nullptr, nullptr);

  // check whether derivative changes sign
  if (GetSign(derivative_position_error_at_param_min) !=
      GetSign(derivative_position_error_at_param_max)) {
    // use bisection method to find the extremum
    double param_extremum = GetRootABKMethod(dp_func, param_min_max);

    // position error at extremum
    double execution_time_at_param_extremum = 0.0;
    double position_error_at_param_extremum =
        p_func(param_extremum, nullptr, &execution_time_at_param_extremum);

    if (position_error_at_param_extremum < 0.0) {
      if (execution_time_at_param_extremum > execution_time) {
        if (param_extremum >= 0.0) {
          if (param_extremum > param_min_max.min) {
            return {param_min_max.min, param_extremum};
          }
        } else {
          if (param_extremum < param_min_max.max) {
            return {param_extremum, param_min_max.max};
          }
        }
      }
    }
  }
  return param_min_max;
}

constexpr auto UpdateLimitsForPosTri = Curry(UpdateLimitsForPosTriFunc);

Range UpdateLimitsForNegLinPosTri(const bool get_min_root,
                                  const Range& param_min_max,
                                  const MotionState& s) {
  // Epsilon value to check whether the function value of an interval limit
  // equals zero.
  constexpr double kBorderlineEpsilon = 1.0e-4;

  constexpr double kUpperLimitDecrementRatio = 1.0e-2;

  // Max number of iterations of upper limit decrement
  constexpr unsigned int kUpperLimitDecrementMaxIterations = 50;

  auto p_func = GetPositionErrorForNegLinPosTriFunc(s);
  auto dp_func = GetDerivativePositionErrorForNegLinPosTri(s);
  auto ddp_func = GetSecondDerivativePositionErrorForNegLinPosTri(s);

  if (param_min_max.min == param_min_max.max) {
    return param_min_max;
  }

  Range result = param_min_max;

  double position_error_at_param_min = p_func(result.min, nullptr, nullptr);
  double position_error_at_param_max = p_func(result.max, nullptr, nullptr);

  double derivative_position_error_at_param_min =
      dp_func(result.min, nullptr, nullptr);
  double derivative_position_error_at_param_max =
      dp_func(result.max, nullptr, nullptr);

  double second_derivative_position_error_at_param_min =
      ddp_func(result.min, nullptr, nullptr);
  double second_derivative_position_error_at_param_max =
      ddp_func(result.max, nullptr, nullptr);

  // this is to avoid a tricky case that requires third derivative
  double new_param_max = result.max;
  unsigned int upper_limit_decrement_iter = 0;
  bool upper_limit_modified = false;
  while (new_param_max > result.min &&
         (GetSign(second_derivative_position_error_at_param_min) ==
          GetSign(second_derivative_position_error_at_param_max)) &&
         (upper_limit_decrement_iter < kUpperLimitDecrementMaxIterations)) {
    new_param_max -= kUpperLimitDecrementRatio * (result.max - result.min);
    second_derivative_position_error_at_param_max =
        ddp_func(new_param_max, nullptr, nullptr);
    upper_limit_modified = true;
    upper_limit_decrement_iter++;
  }

  if (new_param_max < result.min ||
      upper_limit_decrement_iter >= kUpperLimitDecrementMaxIterations) {
    // something fishy - reset to old values
    new_param_max = result.max;
    second_derivative_position_error_at_param_max =
        ddp_func(new_param_max, nullptr, nullptr);
    upper_limit_modified = false;
  }
  result.max = new_param_max;

  if (upper_limit_modified) {
    position_error_at_param_max = p_func(result.max, nullptr, nullptr);
    derivative_position_error_at_param_max =
        dp_func(result.max, nullptr, nullptr);
  }

  // check whether derivative changes sign
  if (GetSign(derivative_position_error_at_param_min) !=
      GetSign(derivative_position_error_at_param_max)) {
    // use bisection method to find the extremum
    double param_extremum = GetRootABKMethod(dp_func, {result.min, result.max});

    // position error at extremum
    double position_error_at_param_extremum =
        p_func(param_extremum, nullptr, nullptr);

    if (get_min_root) {
      if ((GetSign(position_error_at_param_extremum) !=
           GetSign(position_error_at_param_min)) ||
          (EpsilonEqual(position_error_at_param_max, 0.0,
                        kBorderlineEpsilon))) {
        if (param_extremum > result.min) {
          result.max = param_extremum;
          return result;
        }
      }
    } else {
      if ((GetSign(position_error_at_param_extremum) !=
           GetSign(position_error_at_param_max)) ||
          (EpsilonEqual(position_error_at_param_extremum, 0.0,
                        kBorderlineEpsilon))) {
        if (param_extremum < result.max) {
          result.min = param_extremum;
          return result;
        }
      }
    }
  } else {
    // use bisection method to find the extremum of derivative position error
    double param_dp_extremum = GetRootABKMethod(ddp_func, result);

    // derivative of position error at extremum
    double derivative_position_error_at_param_dp_extremum =
        dp_func(param_dp_extremum, nullptr, nullptr);

    // check whether derivative changes sign
    if (GetSign(derivative_position_error_at_param_min) !=
        GetSign(derivative_position_error_at_param_dp_extremum)) {
      // use bisection method to find the two extrema of position error
      double param_lower_extremum =
          GetRootABKMethod(dp_func, {result.min, param_dp_extremum});
      double param_upper_extremum =
          GetRootABKMethod(dp_func, {param_dp_extremum, result.max});

      // position error at extrema
      double position_error_at_param_lower_extremum =
          p_func(param_lower_extremum, nullptr, nullptr);
      double position_error_at_param_upper_extremum =
          p_func(param_upper_extremum, nullptr, nullptr);

      if (get_min_root) {
        if (((GetSign(position_error_at_param_lower_extremum) !=
              GetSign(position_error_at_param_upper_extremum)) ||
             (EpsilonEqual(position_error_at_param_upper_extremum, 0.0,
                           kBorderlineEpsilon))) &&
            (GetSign(position_error_at_param_min) !=
             GetSign(position_error_at_param_lower_extremum))) {
          if (param_lower_extremum > result.min) {
            result.max = param_lower_extremum;
            return result;
          }
        }
      } else {
        if ((GetSign(position_error_at_param_lower_extremum) !=
             GetSign(position_error_at_param_upper_extremum)) ||
            (EpsilonEqual(position_error_at_param_lower_extremum, 0.0,
                          kBorderlineEpsilon))) {
          if (param_lower_extremum < param_upper_extremum) {
            return {param_lower_extremum, param_upper_extremum};
          }
        }
      }
    }
  }
  return result;
}

// Template Check whether a valid solution exists using the chosen profile
bool IsThereAValidSolutionNegLinPosTri(Range param_min_max,
                                       const MotionState& s) {
  bool invalid_solution_min = false;
  bool invalid_solution_max = false;

  auto p_func = GetPositionErrorForNegLinPosTriFunc(s);
  auto dp_func = GetDerivativePositionErrorForNegLinPosTri(s);

  double p_of_min = p_func(param_min_max.min, &invalid_solution_min, nullptr);

  double p_of_max = p_func(param_min_max.max, &invalid_solution_max, nullptr);

  if (GetSign(p_of_min) != GetSign(p_of_max)) {
    return true;
  }

  double dp_of_min = dp_func(param_min_max.min, nullptr, nullptr);
  double dp_of_max = dp_func(param_min_max.max, nullptr, nullptr);

  // Special case when derivative of min/max is close to zero.
  if (fabs(dp_of_min) < kAbsoluteLimitOffset ||
      fabs(dp_of_max) < kAbsoluteLimitOffset) {
    // shrink the range of values
    param_min_max = ShrinkRange(param_min_max);
    dp_of_min = dp_func(param_min_max.min, nullptr, nullptr);
    dp_of_max = dp_func(param_min_max.max, nullptr, nullptr);
  }

  if (GetSign(dp_of_min) != GetSign(dp_of_max)) {
    double param_extremum = GetRootABKMethod(dp_func, param_min_max);

    bool invalid_solution_extremum = false;
    double p_of_extremum =
        p_func(param_extremum, &invalid_solution_extremum, nullptr);

    if (GetSign(p_of_extremum) != GetSign(p_of_max)) {
      return true;
    } else {
      if (fabs(p_of_extremum) < kValidSolutionEpsilon) {
        return !invalid_solution_extremum;
      }
      return false;
    }
  } else {
    if ((fabs(p_of_min) < kValidSolutionEpsilon) ||
        (fabs(p_of_max) < kValidSolutionEpsilon)) {
      if (invalid_solution_min || invalid_solution_max) {
        return false;
      } else {
        return true;
      }
    } else {
      return false;
    }
  }
}

inline MotionState CalcExecTimePosTriNegTriFunc(
    const double step1b_red_profile_execution_time, const MotionState& s) {
  return CalcExecTimePosTriNegTri(UpdateLimitsForPosTri,
                                  step1b_red_profile_execution_time, s);
}
constexpr auto CalcExecTimePosTriNegTri = Curry(CalcExecTimePosTriNegTriFunc);

inline MotionState CalcExecTimeNegTriPosTriFunc(
    const double step1a_red_profile_execution_time, const MotionState& s) {
  return FlipState(s) |
         CalcExecTimePosTriNegTri(step1a_red_profile_execution_time);
}
constexpr auto CalcExecTimeNegTriPosTri = Curry(CalcExecTimeNegTriPosTriFunc);

std::pair<MotionState, double> CalcExecTimeNegLinPosTriWithRootFunc(
    const bool is_step1b3, const bool get_min_root, const MotionState& s) {
  // limits for param
  auto param_min_max = GetNegLinEndAccelerationLimitsForNegLinPosTri(s);

  if (is_step1b3) {
    param_min_max = UpdateLimitsForNegLinPosTri(get_min_root, param_min_max, s);
  }

  // return execution time
  return CalcExecTimeWithRoot(Profile::kNegLinPosTri,
                              GetPositionErrorForNegLinPosTriFunc(s),
                              GetDerivativePositionErrorForNegLinPosTri(s),
                              get_min_root, param_min_max, s);
}

constexpr auto CalcExecTimeNegLinPosTriWithRoot =
    Curry(CalcExecTimeNegLinPosTriWithRootFunc);

constexpr auto CalcExecTimeNegLinPosTri =
    Curry<const bool, const bool, const MotionState&>(
        [](const bool is_step1b3, const bool get_min_root,
           const MotionState& s) {
          return CalcExecTimeNegLinPosTriWithRootFunc(is_step1b3, get_min_root,
                                                      s)
              .first;
        });

MotionState CalcExecTimeMinNegLinPosTrapNegLinPosTriFunc(
    const double step1a_red_profile_execution_time, const MotionState& s) {
  auto s_neglinpostri = CalcExecTimeNegLinPosTri(/*is_step1b3*/ false,
                                                 /*get_min_root*/ true, s);
  auto s_neglinpostrap = CalcExecTimeNegLinPosTrap(
      /*get_min_root*/ true, /*try_limits_if_needed*/ false,
      step1a_red_profile_execution_time, s);

  // Return neglinpostri if it is the only one to succeed, or it has a valid and
  // lower execution time than neglinpostrap.  Otherwise if it fails, or both
  // fails, or it has an invalid time, return neglinpostrap.
  if (IsSuccess(s_neglinpostri) &&
      (IsFailure(s_neglinpostrap) ||
       ((GetT(s_neglinpostri) - GetT(s) > step1a_red_profile_execution_time) &&
        (GetT(s_neglinpostri) < GetT(s_neglinpostrap))))) {
    return s_neglinpostri;
  }

  // NegLinPosTrap is the solution
  return s_neglinpostrap;
}

constexpr auto CalcExecTimeMinNegLinPosTrapNegLinPosTri =
    Curry(CalcExecTimeMinNegLinPosTrapNegLinPosTriFunc);

// Analyze the NegLinPosTrap profile, Calculate one or (if possible) both
// solutions to check for an inoperative time interval. Is one existent?
std::pair<MotionState, MotionState> CalcDecision69(
    const double step1a_red_profile_execution_time, const MotionState& s) {
  auto fail_result =
      std::make_pair(FinishInfiniteFailure(s), FinishInfiniteFailure(s));

  // get execution time, error and root values
  auto state_for_step_b = CalcExecTimeNegLinPosTrap(
      /*get_min_root*/ false, /*try_limits_if_needed*/ false,
      step1a_red_profile_execution_time, s);
  double step_b_delta_t = GetT(state_for_step_b) - GetT(s);
  if (IsFailure(state_for_step_b) ||
      step_b_delta_t <= step1a_red_profile_execution_time ||
      EpsilonEqual(step_b_delta_t, step1a_red_profile_execution_time,
                   kStep1TimeEpsilon)) {
    return fail_result;
  }

  auto state_for_step_c = CalcExecTimeNegLinPosTrap(
      /*get_min_root*/ true, /*try_limits_if_needed*/ false, step_b_delta_t, s);
  if (IsSuccess(state_for_step_b) &&
      GetT(state_for_step_c) > GetT(state_for_step_b) &&
      !EpsilonEqual(GetT(state_for_step_b), GetT(state_for_step_c),
                    kStep1TimeEpsilon)) {
    return std::make_pair(state_for_step_b, state_for_step_c);
  }

  return fail_result;
}

// Analyze the NegLinPosTrap profile, Calculate one or (if possible) both
// solutions to check for an inoperative time interval and also consider
// NegLinPosTri. Is one existent?
std::pair<MotionState, MotionState> CalcDecision71(
    const double step1a_red_profile_execution_time, const MotionState& s) {
  auto fail_result = std::make_pair(SetFailure(s) | SetT(kInfinity),
                                    SetFailure(s) | SetT(kInfinity));

  // get execution time, error and root values
  auto state_for_step_b = CalcExecTimeNegLinPosTrap(
      /*get_min_root*/ false, /*try_limits_if_needed*/ false,
      step1a_red_profile_execution_time, s);
  if (IsFailure(state_for_step_b)) {
    return fail_result;
  }

  auto state_for_step_c_postrap = CalcExecTimeNegLinPosTrap(
      /*get_min_root*/ true, /*try_limits_if_needed*/ false,
      GetT(state_for_step_b) - GetT(s), s);
  if (IsFailure(state_for_step_c_postrap)) {
    return fail_result;
  }

  // check if the start operating zone solution was found
  if (GetT(state_for_step_c_postrap) > GetT(state_for_step_b) &&
      !EpsilonEqual(GetT(state_for_step_b), GetT(state_for_step_c_postrap),
                    kStep1TimeEpsilon)) {
    return std::make_pair(state_for_step_b, state_for_step_c_postrap);
  }

  auto state_for_step_c_postri = CalcExecTimeNegLinPosTri(
      /*is_step1b3*/ true, /*get_min_root*/ true, s);

  // check if the start time solution was found
  if (IsFailure(state_for_step_c_postri) ||
      GetT(state_for_step_c_postri) <= GetT(state_for_step_b) ||
      EpsilonEqual(GetT(state_for_step_b), GetT(state_for_step_c_postri),
                   kStep1TimeEpsilon)) {
    return fail_result;
  }

  return std::make_pair(state_for_step_b, state_for_step_c_postri);
}

// Analyze the NegLinPosTri profile, Calculate one or (if possible) both
// solutions to check for an inoperative time interval. Is one existent?
std::pair<MotionState, MotionState> CalcDecision72(
    const double step1a_red_profile_execution_time, const MotionState& s) {
  auto fail_result = std::make_pair(SetFailure(s) | SetT(kInfinity),
                                    SetFailure(s) | SetT(kInfinity));

  auto param_min_max = UpdateLimitsForNegLinPosTri(
      false, GetNegLinEndAccelerationLimitsForNegLinPosTri(s), s);

  if (!IsThereAValidSolutionNegLinPosTri(param_min_max, s)) {
    return fail_result;
  }

  // get execution time, error and root values
  auto state_for_step_b = CalcExecTimeNegLinPosTri(
      /*is_step1b3*/ true, /*get_min_root*/ false, s);
  double execution_time_step1b = GetT(state_for_step_b) - GetT(s);

  if (IsFailure(state_for_step_b)) {
    return fail_result;
  }

  auto state_for_step_c = CalcExecTimeNegLinPosTri(
      /*is_step1b3*/ true, /*get_min_root*/ true, s);
  if (IsFailure(state_for_step_c)) {
    return fail_result;
  }

  // this is to handle a super ugly case where 1A solution = 1B solution
  if (execution_time_step1b < step1a_red_profile_execution_time ||
      EpsilonEqual(GetT(state_for_step_b), GetT(state_for_step_c),
                   kStep1TimeEpsilon)) {
    auto state_for_step_b_fixed =
        IncrementT(step1a_red_profile_execution_time, s) |
        SetProfile(GetStep1AProfile(s)) | SetSuccess;
    return std::make_pair(state_for_step_b_fixed, state_for_step_c);
  }

  return std::make_pair(state_for_step_b, state_for_step_c);
}

// Analyze the NegLinPosTri profile, Calculate one or (if possible) both
// solutions to check for an inoperative time interval. Is one existent?
std::pair<MotionState, MotionState> CalcDecision83(
    const double step1a_red_profile_execution_time, const MotionState& s) {
  auto fail_result = std::make_pair(SetFailure(s) | SetT(kInfinity),
                                    SetFailure(s) | SetT(kInfinity));

  auto param_min_max = UpdateLimitsForNegLinPosTri(
      false, GetNegLinEndAccelerationLimitsForNegLinPosTri(s), s);

  if (!IsThereAValidSolutionNegLinPosTri(param_min_max, s)) {
    return fail_result;
  }
  // get execution time, error and root values

  auto [state_for_step_b, step_b_root] = CalcExecTimeNegLinPosTriWithRoot(
      /*is_step1b3*/ true, /*get_min_root*/ false, s);
  double s1b_dt = GetT(state_for_step_b) - GetT(s);
  if (IsFailure(state_for_step_b) ||
      s1b_dt < step1a_red_profile_execution_time ||
      !EpsilonEqual(s1b_dt, step1a_red_profile_execution_time,
                    kStep1TimeEpsilon)) {
    return fail_result;
  }

  auto [state_for_step_c, step_c_root] = CalcExecTimeNegLinPosTriWithRoot(
      /*is_step1b3*/ true, /*get_min_root*/ true, s);
  if (IsFailure(state_for_step_c) ||
      GetT(state_for_step_c) <= GetT(state_for_step_b) ||
      EpsilonEqual(GetT(state_for_step_b), GetT(state_for_step_c),
                   kStep1TimeEpsilon)) {
    return fail_result;
  }

  // The following code ensures the consistency and validity of the derived
  // solutions

  // Find the mean root value of 1B and 1C and its position error
  // If the position error at mean is insignificant, then the 1B and 1C
  // solutions are identical, which implies there is no inoperative time
  // interval and they should be infinity
  double mean_root_value = 0.5 * (step_b_root + step_c_root);
  double position_error_at_mean_root_value =
      GetPositionErrorForNegLinPosTriFunc(s, mean_root_value, nullptr, nullptr);

  // Find the position error at half the root value of 1B
  // If a valid inoperative time interval exists then this error must be
  // significantly larger than zero and positive
  double half_root_value_step1b = 0.5 * step_b_root;
  double position_error_at_half_root_value_step1b =
      GetPositionErrorForNegLinPosTriFunc(s, half_root_value_step1b, nullptr,
                                          nullptr);

  if ((position_error_at_mean_root_value <
       kAbsoluteStep1B1Epsilon +
           kRelativeStep1B1Epsilon * fabs(GetP(s) - GetPTrgt(s))) &&
      (position_error_at_half_root_value_step1b <
       kAbsoluteStep1B1Epsilon +
           kRelativeStep1B1Epsilon * fabs(GetP(s) - GetPTrgt(s)))) {
    return fail_result;
  }

  // Find the position error at half the root value of 1C
  // If 1C solution is valid then this error must be less than zero
  double half_root_value_step1c = 0.5 * step_c_root;
  double position_error_at_half_root_value_step1c =
      GetPositionErrorForNegLinPosTriFunc(s, half_root_value_step1c, nullptr,
                                          nullptr);

  if (position_error_at_half_root_value_step1c >= 0.0) {
    return std::make_pair(state_for_step_b, SetT(kInfinity, state_for_step_c));
  }

  return std::make_pair(state_for_step_b, state_for_step_c);
}

template <typename FalseFunc>
MotionState TryApplyCalcExecTimeBC(
    int decision_num, const std::pair<MotionState, MotionState>& sb_and_sc,
    const FalseFunc& ffunc, const MotionState& s) {
  auto [state_for_step_b, state_for_step_c] = sb_and_sc;
  if (IsSuccess(state_for_step_b)) {
    // update Step1C result as well since we have found it
    // already
    auto& step1c = GetDOFOutput(s).step1c;
    WriteStateToSubStepOutput(state_for_step_c, step1c);
    GetProfileTrace(s) = step1c.applied_profile_trace;
    GetDOFOutput(s).inoperative_end_execution_time = GetT(state_for_step_c);

    return state_for_step_b | AppendDecisionTrace(decision_num, true);
  }
  return AppendDecisionTrace(decision_num, false, s) | ffunc;
}

double GetCalcExecTimeExecutionTime(const MotionState& s) {
  return GetDOFOutput(s).min_execution_time - GetT(s);
}

// The decision tree.
// This is done as a struct to allow the decisions to be listed in order
// starting from the first decision.  Bundling the decisions into a struct
// allows us to refer to later decisions before they are declared, whereas a
// namespace wouldn't allow that.
struct Decisions {
  inline static MotionState D1(const MotionState& s) {
    return MakeDecision(1, CmpAGteZero, &D2, FlipState | &D2, s);
  }

  inline static MotionState D2(const MotionState& s) {
    return MakeDecision(2, CmpALteAMax, &D3,
                        CalcExecTimeNegLinAdowntoAMax | &D3, s);
  }

  inline static MotionState D3(const MotionState& s) {
    return MakeDecision(3, ADownToZero | CmpVLteVMax, &D4,
                        CalcExecTimeNegLinAdowntoZero | FlipState |
                            IfElse(GetIntoBoundariesFast(s), &D5, &D87),
                        s);
  }

  inline static MotionState D4(const MotionState& s) {
    return MakeDecision(4, CmpVGteVMin, &D9,
                        IfElse(GetIntoBoundariesFast(s), &D5, &D87), s);
  }

  inline static MotionState D5(const MotionState& s) {
    return MakeDecision(5, AUpToAMax | CmpVGtVMin, &D6, &D7, s);
  }

  inline static MotionState D6(const MotionState& s) {
    return MakeDecision(6, VUpToVMinPosLin | ADownToZero | CmpVLteVMax,
                        CalcExecTimePosLinVuptoVMin | &D9,
                        CalcExecTimePosLinNegLinVuptoVMin | &D9, s);
  }

  inline static MotionState D7(const MotionState& s) {
    return MakeDecision(
        7, SetV(GetVMin(s), s) | SetA(GetAMax(s)) | AMaxDownToZero | CmpVGtVMax,
        &D8, CalcExecTimePosLinHldVuptoVMin | &D9, s);
  }

  inline static MotionState D8(const MotionState& s) {
    return MakeDecision(8, AUpToAMaxDownToZero | CmpVGtVMax,
                        CalcExecTimePosLinNegLinVuptoVMin | &D9,
                        CalcExecTimePosLinHldNegLinVuptoVMin | &D9, s);
  }

  inline static MotionState D9(const MotionState& s) {
    return MakeDecision(9, ADownToZero | CmpVSignEqVTrgt, &D78, &D64, s);
  }

  inline static MotionState D11(const MotionState& s) {
    return MakeDecision(11, CmpAGteZero, &D12, &D48, s);
  }

  inline static MotionState D12(const MotionState& s) {
    return MakeDecision(12, ADownToZero | CmpVLteVTrgt, &D13, &D37, s);
  }

  inline static MotionState D13(const MotionState& s) {
    return MakeDecision(13, AUpToAMaxDownToZero | CmpVLteVTrgt, &D14, &D36, s);
  }

  inline static MotionState D14(const MotionState& s) {
    return MakeDecision(14, VUpToVTrgtPosTrap | CmpPLteEpsilonPTrgt, &D15,
                        FinishInfinite, s);
  }

  inline static MotionState D15(const MotionState& s) {
    return MakeDecision(15, ADownToZero | AUpToAMaxDownToZero | CmpVLteVTrgt,
                        &D16, &D24, s);
  }

  inline static MotionState D16(const MotionState& s) {
    return MakeDecision(
        16, ADownToZero | VUpToVTrgtPosTrap | CmpPLtePTrgt,
        CalcExecTimeNegLinAdowntoZero | &D17,
        CalcExecTimeNegLinPosTrap(/*get_min_root*/ true,
                                  /*try_limits_if_needed*/ false,
                                  GetCalcExecTimeExecutionTime(s)) |
            IfErrorFinishInfinite,
        s);
  }

  inline static MotionState D17(const MotionState& s) {
    return MakeDecision(17, ADownToAMinUpToZero | CmpVGteZero, &D19, &D23, s);
  }

  inline static MotionState D19(const MotionState& s) {
    return TryApplyProfileFunction(
        19,
        CalcExecTimeNegTriPosTrap(/*get_min_root*/ true,
                                  /*try_limits_if_needed*/ false,
                                  GetCalcExecTimeExecutionTime(s)),
        &D21, s);
  }

  inline static MotionState D21(const MotionState& s) {
    // NegTrapPosTrap is the flipped version of PosTrapNegTrap
    return TryApplyProfileFunction(
        19,
        FlipState | CalcExecTimePosTrapNegTrap(/*get_min_root*/ true,
                                               /*try_limits_if_needed*/ false,
                                               GetCalcExecTimeExecutionTime(s)),
        FinishInfinite, s);
  }

  inline static MotionState D23(const MotionState& s) {
    return TryApplyProfileFunction(
        19,
        CalcExecTimeNegTriPosTrap(/*get_min_root*/ true,
                                  /*try_limits_if_needed*/ false,
                                  GetCalcExecTimeExecutionTime(s)),
        FinishInfinite, s);
  }

  inline static MotionState D24(const MotionState& s) {
    return MakeDecision(24,
                        ADownToZero | VUpToVTrgtPosTri | CmpPLteEpsilonPTrgt,
                        CalcExecTimeNegLinAdowntoZero | &D25, &D35, s);
  }

  inline static MotionState D25(const MotionState& s) {
    return MakeDecision(25, ADownToAMinUpToZero | CmpVGteZero, &D75, &D29, s);
  }

  inline static MotionState D26(const MotionState& s) {
    // if a->amin->0 and then v->vtrgt (PosTrap), is v>=0?
    bool decision =
        FlipState(s) | AUpToAMaxDownToZero | VDownToVTrgtNegTrap | CmpPLtePTrgt;
    return MakeDecision(26, decision, &D27, &D28, s);
  }

  inline static MotionState D27(const MotionState& s) {
    // if a->amin->0 and then a->amax->0 so that v=vtrgt, is p>=ptrgt?
    // Input is flipped, so this becomes a->apeak(<amax)->0->amin->0
    return MakeDecision(
        27, VToVTrgtPosTriAndThen(ADownToAMinUpToZero) | CmpPLtePTrgt,
        CalcExecTimePosTriNegTri(GetCalcExecTimeExecutionTime(s)),
        CalcExecTimePosTriNegTrap(/*get_min_root*/ true,
                                  /*try_limits_if_needed*/ false,
                                  GetCalcExecTimeExecutionTime(s)) |
            IfErrorFinishInfinite,
        FlipState(s));
  }

  inline static MotionState D28(const MotionState& s) {
    return TryApplyProfileFunction(
        28, CalcExecTimeNegTriPosTri(GetCalcExecTimeExecutionTime(s)), &D19, s);
  }

  inline static MotionState D29(const MotionState& s) {
    // if a->-apeak->0, then set vx=v, and then a->amax->0 so that v=vtrgt, is
    // vx>=0?
    return MakeDecision(
        29,
        (GetVTrgt(s) - GetV(SetV(0, s) | SetA(0) | AUpToAMaxDownToZero)) >= 0,
        &D30, &D34, s);
  }

  inline static MotionState D30(const MotionState& s) {
    return MakeDecision(30,
                        VDownToZeroNegTri | VUpToVTrgtPosTrap | CmpPGtePTrgt,
                        &D27, &D31, s);
  }

  inline static MotionState D31(const MotionState& s) {
    return TryApplyProfileFunction(
        31, CalcExecTimeNegTriPosTri(GetCalcExecTimeExecutionTime(s)), &D32, s);
  }

  inline static MotionState D32(const MotionState& s) {
    return MakeDecision(
        19,
        CalcExecTimeNegTriPosTrap(/*get_min_root*/ true,
                                  /*try_limits_if_needed*/ false,
                                  GetCalcExecTimeExecutionTime(s)) |
            IsFailure,
        FinishInfinite, &D27, s);
  }

  inline static MotionState D34(const MotionState& s) {
    return TryApplyProfileFunction(
        34, CalcExecTimeNegTriPosTri(GetCalcExecTimeExecutionTime(s)),
        FinishInfinite, s);
  }

  inline static MotionState D35(const MotionState& s) {
    return CalcExecTimeMinNegLinPosTrapNegLinPosTri(
               GetCalcExecTimeExecutionTime(s), s) |
           IfErrorFinishInfinite;
  }

  inline static MotionState D36(const MotionState& s) {
    return MakeDecision(36, VUpToVTrgtPosTri | CmpPLteEpsilonPTrgt, &D24,
                        FinishInfinite, s);
  }

  inline static MotionState D37(const MotionState& s) {
    return MakeDecision(37, ADownToZero | ADownToAMinUpToZero | CmpVLteVTrgt,
                        &D38, &D47, s);
  }

  inline static MotionState D38(const MotionState& s) {
    return MakeDecision(
        38, ADownToZero | VDownToVTrgtNegTri | CmpPLteEpsilonPTrgt,
        CalcExecTimeNegLinAdowntoZero | &D39, FinishInfinite, s);
  }

  inline static MotionState D39(const MotionState& s) {
    return MakeDecision(39, ADownToAMinUpToZero | CmpVGteZero, &D75, &D76, s);
  }

  inline static MotionState D41(const MotionState& s) {
    return TryApplyProfileFunction(
        41, CalcExecTimeNegTriPosTri(GetCalcExecTimeExecutionTime(s)), &D42, s);
  }

  inline static MotionState D42(const MotionState& s) {
    // if a->amin->hold->0, then set vx=v, and then a->amax->0 so that v=vtrgt,
    // is vx>=0?
    return MakeDecision(
        42,
        (GetVTrgt(s) - GetV(SetV(0, s) | SetA(0) | AUpToAMaxDownToZero) >= 0),
        &D44, &D46, s);
  }

  inline static MotionState D44(const MotionState& s) {
    return TryApplyProfileFunction(
        44,
        CalcExecTimeNegTrapPosTri(/*get_min_root*/ false,
                                  /*try_limits_if_needed*/ false,
                                  GetCalcExecTimeExecutionTime(s)),
        &D21, s);
  }

  inline static MotionState D46(const MotionState& s) {
    return TryApplyProfileFunction(
        44,
        CalcExecTimeNegTrapPosTri(/*get_min_root*/ false,
                                  /*try_limits_if_needed*/ false,
                                  GetCalcExecTimeExecutionTime(s)),
        FinishInfinite, s);
  }

  inline static MotionState D47(const MotionState& s) {
    return MakeDecision(
        47, ADownToZero | VDownToVTrgtNegTrap | CmpPLteEpsilonPTrgt,
        CalcExecTimeNegLinAdowntoZero | &D42, FinishInfinite, s);
  }

  inline static MotionState D48(const MotionState& s) {
    return MakeDecision(48, AUpToZero | CmpVGteVTrgt, &D49, &D57, s);
  }

  inline static MotionState D49(const MotionState& s) {
    return MakeDecision(49, ADownToAMinUpToZero | CmpVGteVTrgt, &D56, &D50, s);
  }

  inline static MotionState D50(const MotionState& s) {
    return MakeDecision(50, VDownToVTrgtNegTri | CmpPGteEpsilonPTrgt,
                        FinishInfinite, &D53, s);
  }

  inline static MotionState D51(const MotionState& s) {
    return MakeDecision(51, AUpToZero | ADownToAMinUpToZero | CmpVGteVTrgt,
                        &D52, &D53, s);
  }

  inline static MotionState D52(const MotionState& s) {
    return MakeDecision(
        52, AUpToZero | VDownToVTrgtNegTrap | CmpPGtePTrgt,
        FlipState |
            CalcExecTimeNegLinPosTrap(
                /*get_min_root*/ true, /*try_limits_if_needed*/ false,
                GetCalcExecTimeExecutionTime(s)) |
            IfErrorFinishInfinite,
        &D42, s);
  }

  inline static MotionState D53(const MotionState& s) {
    return MakeDecision(53, AUpToZero | VDownToVTrgtNegTri | CmpPGtePTrgt,
                        FlipState | &D54, &D55, s);
  }

  inline static MotionState D54(const MotionState& s) {
    return CalcExecTimeMinNegLinPosTrapNegLinPosTri(
               GetCalcExecTimeExecutionTime(s), s) |
           IfErrorFinishInfinite;
  }

  inline static MotionState D55(const MotionState& s) {
    return MakeDecision(55, ADownToAMinUpToZero | CmpVGteVTrgt, &D42, &D39, s);
  }

  inline static MotionState D56(const MotionState& s) {
    return MakeDecision(56, VDownToVTrgtNegTrap | CmpPGteEpsilonPTrgt,
                        FinishInfinite, &D51, s);
  }

  inline static MotionState D57(const MotionState& s) {
    return MakeDecision(57, AUpToZero | AUpToAMaxDownToZero | CmpVGteVTrgt,
                        &D59, &D58, s);
  }

  inline static MotionState D58(const MotionState& s) {
    return MakeDecision(58, AUpToZero | VUpToVTrgtPosTrap | CmpPGteEpsilonPTrgt,
                        FinishInfinite, &D17, s);
  }

  inline static MotionState D59(const MotionState& s) {
    return MakeDecision(59, AUpToZero | VUpToVTrgtPosTri | CmpPGteEpsilonPTrgt,
                        FinishInfinite, &D60, s);
  }

  inline static MotionState D60(const MotionState& s) {
    return MakeDecision(
        60, ADownToAMinUpToZero | AUpToAMaxDownToZero | CmpVGteVTrgt, &D39,
        &D61, s);
  }

  inline static MotionState D61(const MotionState& s) {
    return MakeDecision(
        61,
        (GetVTrgt(s) - GetV(SetV(0, s) | SetA(0) | AUpToAMaxDownToZero)) >= 0,
        &D63, &D34, s);
  }

  inline static MotionState D63(const MotionState& s) {
    return TryApplyProfileFunction(
        63, CalcExecTimeNegTriPosTri(GetCalcExecTimeExecutionTime(s)), &D17, s);
  }

  inline static MotionState D64(const MotionState& s) {
    return MakeDecision(64, CmpVTrgtLtZero, FinishInfinite, &D65, s);
  }

  inline static MotionState D65(const MotionState& s) {
    // Was the result of Step1A one of the following profiles: PosTriNegTri,
    // PosTrapNegTri, PosTriZeroNegTri, PosTrapZeroNegTri, PosTriNegTrap,
    // PosTrapNegTrap, PosTriZeroNegTrap or PosTrapZeroNegTrap?
    const Profile step1a_profile = GetStep1AProfile(s);
    bool decision = (((step1a_profile == Profile::kPosTriNegTri ||
                       step1a_profile == Profile::kPosTrapNegTri ||
                       step1a_profile == Profile::kPosTriZeroNegTri ||
                       step1a_profile == Profile::kPosTrapZeroNegTri ||
                       step1a_profile == Profile::kPosTriNegTrap ||
                       step1a_profile == Profile::kPosTrapNegTrap ||
                       step1a_profile == Profile::kPosTriZeroNegTrap ||
                       step1a_profile == Profile::kPosTrapZeroNegTrap) &&
                      !GetStep1AAppliedProfileFlipped(s)) ||
                     step1a_profile == Profile::kNegLinPosTri ||
                     step1a_profile == Profile::kNegLinPosTrap);
    return MakeDecision(65, decision, &D66, FinishInfinite, s);
  }

  inline static MotionState D66(const MotionState& s) {
    return MakeDecision(66, AUpToAMaxDownToZero | CmpVLteVTrgt, &D67, &D74, s);
  }

  inline static MotionState D67(const MotionState& s) {
    return MakeDecision(67, ADownToZero | AUpToAMaxDownToZero | CmpVLteVTrgt,
                        &D68, &D70, s);
  }

  inline static MotionState D68(const MotionState& s) {
    return MakeDecision(
        68, ADownToZero | VUpToVTrgtPosTrap | CmpPLtePTrgt, &D69,
        CalcExecTimeNegLinPosTrap(/*get_min_root*/ true,
                                  /*try_limits_if_needed*/ false,
                                  GetCalcExecTimeExecutionTime(s)) |
            IfErrorFinishInfinite,
        s);
  }

  inline static MotionState D69(const MotionState& s) {
    return TryApplyCalcExecTimeBC(
        69, CalcDecision69(GetCalcExecTimeExecutionTime(s), s), FinishInfinite,
        s);
  }

  inline static MotionState D70(const MotionState& s) {
    return MakeDecision(70, ADownToZero | VUpToVTrgtPosTri | CmpPLtePTrgt, &D71,
                        &D73, s);
  }

  inline static MotionState D71(const MotionState& s) {
    return TryApplyCalcExecTimeBC(
        71, CalcDecision71(GetCalcExecTimeExecutionTime(s), s), &D72, s);
  }

  inline static MotionState D72(const MotionState& s) {
    return TryApplyCalcExecTimeBC(
        72, CalcDecision72(GetCalcExecTimeExecutionTime(s), s), FinishInfinite,
        s);
  }

  inline static MotionState D73(const MotionState& s) {
    return CalcExecTimeMinNegLinPosTrapNegLinPosTri(
               GetCalcExecTimeExecutionTime(s), s) |
           IfErrorFinishInfinite;
  }

  inline static MotionState D74(const MotionState& s) {
    return MakeDecision(74, ADownToZero | VUpToVTrgtPosTri | CmpPLtePTrgt, &D72,
                        CalcExecTimeNegLinPosTri(/*is_step_b3*/ false,
                                                 /*get_min_root*/ true),
                        s);
  }

  inline static MotionState D75(const MotionState& s) {
    return MakeDecision(
        75, ADownToAMinUpToZero | AUpToAMaxDownToZero | CmpVGteVTrgt, &D41,
        &D26, s);
  }

  inline static MotionState D76(const MotionState& s) {
    // if v->0 (NegTri) and then a->amax->0, is v>=vtrgt?
    // We can just start at v=a=0, the end of negtri
    return MakeDecision(
        76, SetV(0, s) | SetA(0) | AUpToAMaxDownToZero | CmpVGteVTrgt, &D34,
        &D77, s);
  }

  inline static MotionState D77(const MotionState& s) {
    return TryApplyProfileFunction(
        77, CalcExecTimeNegTriPosTri(GetCalcExecTimeExecutionTime(s)), &D23, s);
  }

  inline static MotionState D78(const MotionState& s) {
    // Was the result of Step1A NegLinPosTrap or NegLinPosTri?
    return MakeDecision(78,
                        ((GetStep1AProfile(s) == Profile::kNegLinPosTrap) ||
                         (GetStep1AProfile(s) == Profile::kNegLinPosTri)),
                        &D80, &D79, s);
  }

  inline static MotionState D79(const MotionState& s) {
    return MakeDecision(79, CmpVTrgtLtZero, FlipState | &D11, &D11, s);
  }

  inline static MotionState D80(const MotionState& s) {
    // Was the result of Step1A NegLinPosTrap?
    return MakeDecision(80, (GetStep1AProfile(s) == Profile::kNegLinPosTrap),
                        &D81, &D83, s);
  }

  inline static MotionState D81(const MotionState& s) {
    return TryApplyCalcExecTimeBC(
        81, CalcDecision69(GetCalcExecTimeExecutionTime(s), s), &D84, s);
  }

  inline static MotionState D82(const MotionState& s) {
    return TryApplyCalcExecTimeBC(
        82, CalcDecision71(GetCalcExecTimeExecutionTime(s), s), &D83, s);
  }

  inline static MotionState D83(const MotionState& s) {
    return TryApplyCalcExecTimeBC(
        83, CalcDecision83(GetCalcExecTimeExecutionTime(s), s), &D85, s);
  }

  inline static MotionState D84(const MotionState& s) {
    return MakeDecision(84, ADownToZero | AUpToAMaxDownToZero | CmpVLteVTrgt,
                        &D85, &D82, s);
  }

  inline static MotionState D85(const MotionState& s) {
    return MakeDecision(85, CmpVTrgtLtZero, FlipState | &D39, &D86, s);
  }

  inline static MotionState D86(const MotionState& s) {
    return MakeDecision(86, ADownToZero | AUpToAMaxDownToZero | CmpVLteVTrgt,
                        CalcExecTimeNegLinAdowntoZero | &D17,
                        CalcExecTimeNegLinAdowntoZero | &D25, s);
  }

  inline static MotionState D87(const MotionState& s) {
    return MakeDecision(87, AUpToAMaxDownToZero | CmpVLteVMin,
                        CalcExecTimePosTrapVuptoVMin | &D9,
                        CalcExecTimePosTriVuptoVMin | &D9, s);
  }
};

}  // namespace

bool CalcInoperativeZoneStart(const PositionInputs& inputs,
                              const PositionFlags& flags,
                              const MaxDOFFixedVector<bool>& step1a_flipped,
                              PositionOutputs& outputs) {
  for (auto [dof_input, dof_output] :
       Zip(inputs.GetDOFs(), outputs.GetDOFs())) {
    // Check if the dof is available for use and if the constraints are valid.
    if (!dof_input.selected) {
      dof_output.step1b.result = true;
      continue;
    }

    MotionState s = RunPositionStep1DecisionTree(
        Decisions::D1, dof_input, flags, dof_output, dof_output.step1b,
        step1a_flipped[dof_input.index]);

    // update the profile error element
    if (IsFailure(s)) {
      // an error occurred so quit deciding for the rest of the dofs
      // and notify failure
      return false;
    }

    dof_output.inoperative_begin_execution_time = GetT(s);
    WriteStateToSubStepOutput(s, dof_output.step1b);
  }

  return true;
}

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic
