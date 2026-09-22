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

#ifndef INTRINSIC_ICON_REFLEXXES_INTERNAL_POSITION_CALC_EXEC_TIME_FUNCS_H_
#define INTRINSIC_ICON_REFLEXXES_INTERNAL_POSITION_CALC_EXEC_TIME_FUNCS_H_

#include <functional>
#include <utility>

#include "absl/strings/str_format.h"
#include "intrinsic/icon/reflexxes/inputs.h"
#include "intrinsic/icon/reflexxes/internal/intermediate_profile_funcs.h"
#include "intrinsic/icon/reflexxes/internal/motion_state.h"
#include "intrinsic/icon/reflexxes/internal/motion_state_transforms.h"
#include "intrinsic/icon/reflexxes/internal/position_calc_min_execution_time_base.h"
#include "intrinsic/icon/reflexxes/outputs.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

// The (integer) number of maximal loops for the
// Anderson-Bjoerck-King method within the red profiles.
constexpr unsigned int kNumberOfLoopsABKRed = 60;

constexpr auto CalcExecTimePosTrapZeroNegTrap =
    CalcExecTimePosTrapVuptoVMax | HoldToPTrgtAndThen(VDownToVTrgtNegTrap) |
    SetProfile(Profile::kPosTrapZeroNegTrap);

constexpr auto CalcExecTimePosTriZeroNegTrap =
    CalcExecTimePosTriVuptoVMax | HoldToPTrgtAndThen(VDownToVTrgtNegTrap) |
    SetProfile(Profile::kPosTriZeroNegTrap);

constexpr auto CalcExecTimePosTrapZeroNegTri =
    CalcExecTimePosTrapVuptoVMax | HoldToPTrgtAndThen(VDownToVTrgtNegTri) |
    SetProfile(Profile::kPosTrapZeroNegTri);

constexpr auto CalcExecTimePosTriZeroNegTri =
    CalcExecTimePosTriVuptoVMax | HoldToPTrgtAndThen(VDownToVTrgtNegTri) |
    SetProfile(Profile::kPosTriZeroNegTri);

MotionState CalcExecTimePosTrapNegTrapFunc(bool get_min_root,
                                           bool try_limits_if_needed,
                                           double min_execution_time,
                                           const MotionState& s);
constexpr auto CalcExecTimePosTrapNegTrap =
    Curry(CalcExecTimePosTrapNegTrapFunc);

MotionState CalcExecTimePosTrapNegTriFunc(bool get_min_root,
                                          bool try_limits_if_needed,
                                          double min_execution_time,
                                          const MotionState& s);
constexpr auto CalcExecTimePosTrapNegTri = Curry(CalcExecTimePosTrapNegTriFunc);

MotionState CalcExecTimePosTriNegTrapFunc(bool get_min_root,
                                          bool try_limits_if_needed,
                                          double min_execution_time,
                                          const MotionState& s);
constexpr auto CalcExecTimePosTriNegTrap = Curry(CalcExecTimePosTriNegTrapFunc);

MotionState CalcExecTimeNegLinPosTrapFunc(bool get_min_root,
                                          bool try_limits_if_needed,
                                          double min_execution_time,
                                          const MotionState& s);
constexpr auto CalcExecTimeNegLinPosTrap = Curry(CalcExecTimeNegLinPosTrapFunc);

constexpr auto CalcExecTimeNegTriPosTrap =
    Curry<bool, bool, double, MotionState>(
        [](const bool get_min_root, const bool try_limits_if_needed,
           const double step1a_red_profile_execution_time,
           const MotionState& s) {
          // NegTriPosTrap is the flipped version of PosTriNegTrap
          return FlipState(s) |
                 CalcExecTimePosTriNegTrap(get_min_root, try_limits_if_needed,
                                           step1a_red_profile_execution_time);
        });

constexpr auto CalcExecTimeNegTrapPosTrap = Curry<double, MotionState>(
    [](const double step1a_red_profile_execution_time, const MotionState& s) {
      // NegTrapPosTrap is the flipped version of PosTrapNegTrap
      return FlipState(s) | CalcExecTimePosTrapNegTrap(
                                true, false, step1a_red_profile_execution_time);
    });

constexpr auto CalcExecTimeNegTrapPosTri =
    Curry<bool, bool, double, MotionState>(
        [](const bool get_min_root, const bool try_limits_if_needed,
           const double step1a_red_profile_execution_time,
           const MotionState& s) {
          // NegTrapPosTri is the flipped version of PosTrapNegTri
          return FlipState(s) |
                 CalcExecTimePosTrapNegTri(false, false,
                                           step1a_red_profile_execution_time);
        });

// These are a set of limit functions exposes for substeps building custom
// versions of the Red profiles.
Range GetPosTriPeakAccelerationLimitsForPosTriNegTri(const MotionState& s);
Range GetNegLinEndAccelerationLimitsForNegLinPosTri(const MotionState& s);

// Position Error Function for PosTriNegTri
double GetPositionErrorForPosTriNegTriFunc(const MotionState& s,
                                           double peak_acceleration_postri,
                                           bool* invalid_solution,
                                           double* execution_time);

constexpr auto GetPositionErrorForPosTriNegTri =
    Curry(GetPositionErrorForPosTriNegTriFunc);

// Derivative Position Error Function for PosTriNegTri
double GetDerivativePositionErrorForPosTriNegTriFunc(
    const MotionState& s, double peak_acceleration_postri,
    bool* invalid_solution, double* execution_time);

constexpr auto GetDerivativePositionErrorForPosTriNegTri =
    Curry(GetDerivativePositionErrorForPosTriNegTriFunc);

// Position Error Function for NegLinPosTri
double GetPositionErrorForNegLinPosTriFuncFunc(const MotionState& s,
                                               double end_acceleration_neglin,
                                               bool* invalid_solution,
                                               double* execution_time);
constexpr auto GetPositionErrorForNegLinPosTriFunc =
    Curry(GetPositionErrorForNegLinPosTriFuncFunc);

// Derivative of Position Error Function for NegLinPosTrap
double GetDerivativePositionErrorForNegLinPosTriFunc(
    const MotionState& s, double end_acceleration_neglin,
    bool* invalid_solution, double* execution_time);
constexpr auto GetDerivativePositionErrorForNegLinPosTri =
    Curry(GetDerivativePositionErrorForNegLinPosTriFunc);

// Second Derivative of Position Error Function for NegLinPosTrap
double GetSecondDerivativePositionErrorForNegLinPosTriFunc(
    const MotionState& s, double end_acceleration_neglin,
    bool* invalid_solution, double* execution_time);
constexpr auto GetSecondDerivativePositionErrorForNegLinPosTri =
    Curry(GetSecondDerivativePositionErrorForNegLinPosTriFunc);

// Template Anderson-Bjorck-King (ABK) Method to Find Roots
// func is function whose root is to be found, with the form
//    double(double param, bool* invalid_solution, double* execution_time)
// Returns the found root.
template <class F>
double GetRootABKMethod(const F& func, const Range& limits) {
  double param_1 = limits.min;
  double param_2 = limits.max;

  double old_param_1 = 0.0;
  double old_param_2 = 0.0;

  double function_at_param_1 = func(param_1, nullptr, nullptr);

  if (fabs(function_at_param_1) < kABKEpsilon) {
    return param_1;
  }

  double function_at_param_2 = func(param_2, nullptr, nullptr);

  if (fabs(function_at_param_2) < kABKEpsilon) {
    return param_2;
  }

  for (unsigned int i = 0; i < kNumberOfLoopsABKRed; i++) {
    // This part is needed for the case that
    // we can numerically not get closer to the root.

    if ((old_param_1 == param_1) && (old_param_2 == param_2)) {
      return (0.5 * (param_1 + param_2));
    }

    old_param_1 = param_1;
    old_param_2 = param_2;

    if (fabs(param_2 - param_1) <= kABKEpsilon) {
      return (((fabs(function_at_param_2)) < (fabs(function_at_param_1)))
                  ? (param_2)
                  : (param_1));
    }

    // bisection step
    double param_3 = 0.5 * (param_1 + param_2);
    bool invalid_solution = false;
    double function_at_param_3 = func(param_3, &invalid_solution, nullptr);
    if (fabs(function_at_param_3) < kABKEpsilon) {
      return param_3;
    }

    // interval determination step 1
    if ((function_at_param_3 * function_at_param_2) < 0.0) {
      param_1 = param_2;
      param_2 = param_3;
      function_at_param_1 = function_at_param_2;
      function_at_param_2 = function_at_param_3;
    } else {
      param_2 = param_3;
      function_at_param_2 = function_at_param_3;
    }

    // secant step
    double secant_12 =
        (function_at_param_1 - function_at_param_2) / (param_1 - param_2);
    // check for valid secant value
    // ignore the secant step if not valid
    // this SHOULD NOT happen since  the function is monotically increasing or
    // decreasing within the limits
    // but still checking here for completeness
    if (fabs(secant_12) <= kABKFunctionEpsilon) {
      continue;
    }
    param_3 = param_2 - function_at_param_2 / secant_12;
    function_at_param_3 = func(param_3, &invalid_solution, nullptr);
    if (fabs(function_at_param_3) < kABKEpsilon) {
      return param_3;
    }

    // determination of new inclusion interval
    if ((function_at_param_3 * function_at_param_2) < 0.0) {
      param_1 = param_2;
      param_2 = param_3;
      function_at_param_1 = function_at_param_2;
      function_at_param_2 = function_at_param_3;
    } else {
      double g = 1.0 - (function_at_param_3 / function_at_param_2);
      if (g <= 0.0) {
        g = 0.5;
      }
      param_2 = param_3;
      function_at_param_1 *= g;
      function_at_param_2 = function_at_param_3;
    }
  }

  return ((fabs(function_at_param_2) < fabs(function_at_param_1)) ? (param_2)
                                                                  : (param_1));
}

// Provides an additional error tolerance value in order to improve the
// numerical robustness. Takes the time value of the execution time in seconds
// and returns he absolute value of the additionally allowed position error.
double GetAnAdditionalNumericalErrorToleranceBasedOnTheExecutionTime(
    double execution_time_value);

// A base profile functions for a red profile that needs to find a root.
// Takes the position error and derivative position error function of the
// profile. get_min_root is the flag used to decide whether the minimum root
// must be used param_min_max is the limit of the parameter of the chosen
// profile
// Returns the state with the profile applied, and the root value.
template <typename PFunc, typename DPFunc>
std::pair<MotionState, double> CalcExecTimeWithRoot(
    const Profile profile, const PFunc& p_func, const DPFunc& dp_func,
    const bool get_min_root, const Range& param_min_max, const MotionState& s) {
  Range params = param_min_max;
  // original limits
  double original_param_min = params.min;
  double original_param_max = params.max;
  // position error at the original parameter limits
  double position_error_at_original_param_min =
      p_func(params.min, nullptr, nullptr);
  double position_error_at_original_param_max =
      p_func(params.max, nullptr, nullptr);

  // initialize optimal param
  double param_optimal = 0.0;

  // position error at the parameter limits
  double position_error_at_param_min = position_error_at_original_param_min;
  double position_error_at_param_max = position_error_at_original_param_max;

  if (position_error_at_param_min == 0) {
    param_optimal = params.min;
  } else if (position_error_at_param_max == 0) {
    param_optimal = params.max;
  } else {
    // check whether derivative of position error needs to be used
    if (GetSign(position_error_at_param_min) ==
        GetSign(position_error_at_param_max)) {
      // need to compute derivatives
      // derivative of position error at the param_eter limits
      double derivative_position_error_at_param_min =
          dp_func(params.min, nullptr, nullptr);
      double derivative_position_error_at_param_max =
          dp_func(params.max, nullptr, nullptr);

      /////////
      // Special case
      if (fabs(derivative_position_error_at_param_min) < kAbsoluteLimitOffset ||
          fabs(derivative_position_error_at_param_max) < kAbsoluteLimitOffset) {
        // shrink the range of values
        params = ShrinkRange(params);

        position_error_at_param_min = p_func(params.min, nullptr, nullptr);
        position_error_at_param_max = p_func(params.max, nullptr, nullptr);

        derivative_position_error_at_param_min =
            dp_func(params.min, nullptr, nullptr);
        derivative_position_error_at_param_max =
            dp_func(params.max, nullptr, nullptr);
      }
      ////////

      // check whether derivative changes sign
      if (GetSign(derivative_position_error_at_param_min) !=
          GetSign(derivative_position_error_at_param_max)) {
        // use bisection method to find the extremum
        double param_extremum = GetRootABKMethod(dp_func, params);

        // position error at extremum
        double position_error_at_param_extremum =
            p_func(param_extremum, nullptr, nullptr);

        // check whether the sign of the position error at extremum is
        // different from that at min param_eter
        if (GetSign(position_error_at_param_min) !=
            GetSign(position_error_at_param_extremum)) {
          if (get_min_root) {
            // use bisection method to find a root between params.min and
            // param_extremum
            double rootmin =
                GetRootABKMethod(p_func, {params.min, param_extremum});

            bool invalid_solution = false;
            p_func(rootmin, &invalid_solution, nullptr);

            if (invalid_solution) {
              // if the min root returns invalid solution then find max root
              param_optimal =
                  GetRootABKMethod(p_func, {param_extremum, params.max});
            } else {
              param_optimal = rootmin;
            }
          } else {
            // use bisection method to find a root between param_extremum and
            // params.max
            double rootmax =
                GetRootABKMethod(p_func, {param_extremum, params.max});

            bool invalid_solution = false;
            p_func(rootmax, &invalid_solution, nullptr);

            if (invalid_solution) {
              // if the max root returns invalid solution then find min root
              param_optimal =
                  GetRootABKMethod(p_func, {params.min, param_extremum});
            } else {
              param_optimal = rootmax;
            }
          }
        } else {
          // this case occurs due to numerically errors
          if (fabs(position_error_at_original_param_min) <
              fabs(position_error_at_param_min)) {
            params.min = original_param_min;
            position_error_at_param_min = position_error_at_original_param_min;
          }
          if (fabs(position_error_at_original_param_max) <
              fabs(position_error_at_param_max)) {
            params.max = original_param_max;
            position_error_at_param_max = position_error_at_original_param_max;
          }
          if (fabs(position_error_at_param_min) <
              fabs(position_error_at_param_extremum)) {
            if (fabs(position_error_at_param_min) <
                fabs(position_error_at_param_max)) {
              param_optimal = params.min;
            } else {
              param_optimal = params.max;
            }
          } else {
            if (fabs(position_error_at_param_extremum) <
                fabs(position_error_at_param_max)) {
              param_optimal = param_extremum;
            } else {
              param_optimal = params.max;
            }
          }
        }
      } else {
        // this case occurs due to numerically errors
        if (fabs(position_error_at_original_param_min) <
            fabs(position_error_at_param_min)) {
          params.min = original_param_min;
          position_error_at_param_min = position_error_at_original_param_min;
        }
        if (fabs(position_error_at_original_param_max) <
            fabs(position_error_at_param_max)) {
          params.max = original_param_max;
          position_error_at_param_max = position_error_at_original_param_max;
        }
        // return the param corresponding to the lowest position error
        if (fabs(position_error_at_param_min) <
            fabs(position_error_at_param_max)) {
          param_optimal = params.min;
        } else {
          param_optimal = params.max;
        }
      }
    } else {
      // use bisection method to find a root between param limits
      param_optimal = GetRootABKMethod(p_func, params);
    }
  }

  // execution time and positionerror using optimal param
  double execution_time_optimal = 0.0;
  bool invalid_solution = false;
  double position_error_optimal =
      p_func(param_optimal, &invalid_solution, &execution_time_optimal);

  if ((fabs(position_error_optimal) >
       (fabs(kRelStep1PositionErrorTolerance * (GetP(s) - GetPTrgt(s))) +
        kAbsStep1PositionErrorTolerance +
        GetAnAdditionalNumericalErrorToleranceBasedOnTheExecutionTime(
            execution_time_optimal))) ||
      invalid_solution) {
    return std::make_pair(SetFailure(s), 0.0);
  }

  return std::make_pair(
      IncrementT(execution_time_optimal, s) | SetProfile(profile),
      param_optimal);
}

// Version of RedProfile that just returns the state.
template <typename PFunc, typename DPFunc>
MotionState RedProfile(const Profile profile, const PFunc& p_func,
                       const DPFunc& dp_func, const bool get_min_root,
                       const Range& param_min_max, const MotionState& s) {
  return CalcExecTimeWithRoot(profile, p_func, dp_func, get_min_root,
                              param_min_max, s)
      .first;
}

template <typename UpdateLimitsFuncType>
MotionState CalcExecTimePosTriNegTri(
    const UpdateLimitsFuncType& update_limits_func,
    const double min_execution_time, const MotionState& s) {
  auto param_min_max_limit = GetPosTriPeakAccelerationLimitsForPosTriNegTri(s);

  // shrink parameter range for numerical stability
  // update limits for Step1C
  auto s1 =
      RedProfile(Profile::kPosTriNegTri, GetPositionErrorForPosTriNegTri(s),
                 GetDerivativePositionErrorForPosTriNegTri(s), false,
                 update_limits_func(min_execution_time,
                                    ShrinkRange(param_min_max_limit), s),
                 s);

  if (IsSuccess(s1)) {
    return s1;
  }

  // test without shrinking the parameter range
  return RedProfile(
      Profile::kPosTriNegTri, GetPositionErrorForPosTriNegTri(s),
      GetDerivativePositionErrorForPosTriNegTri(s), false,
      update_limits_func(min_execution_time, param_min_max_limit, s), s);
}

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_INTERNAL_POSITION_CALC_EXEC_TIME_FUNCS_H_
