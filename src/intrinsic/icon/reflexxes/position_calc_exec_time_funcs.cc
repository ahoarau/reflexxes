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

#include "intrinsic/icon/reflexxes/internal/position_calc_exec_time_funcs.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include "intrinsic/icon/reflexxes/internal/constants.h"
#include "intrinsic/icon/reflexxes/internal/math.h"
#include "intrinsic/icon/reflexxes/internal/motion_state.h"
#include "intrinsic/icon/reflexxes/internal/motion_state_transforms.h"
#include "intrinsic/icon/reflexxes/internal/polynomial_solvers.h"
#include "intrinsic/icon/reflexxes/internal/position_calc_min_execution_time_base.h"
#include "intrinsic/icon/reflexxes/profile.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

namespace {

// Upper bound of acceptable execution times for position error tolerance.
constexpr double kLargeExecutionTimeLevel1 = 50.0;

// Upper bound of acceptable execution times for position error tolerance
// while iterating to revise the position error tolerance.
constexpr double kLargeExecutionTimeLevel2 = 10.0;

// Ratio by which the execution time is reduced at each iteration
// while revising position error tolerance.
constexpr double kLargeExecutionTimeReductionRatio = 0.1;

// Ratio by which the position error tolerance is increased
// at each iteration during its revision.
constexpr double kAdditionalPositionErrorToleranceIncreaseRatio = 4.0;

// Number of iterations used to revise position error tolerance
// when necessary.
constexpr unsigned int kAdditionalPositionErrorToleranceRevisionIterations = 12;

MotionState FindExecutionTime(
    const Profile profile,
    double (*texe_function)(const MotionState& s, const double root),
    const double min_execution_time, const Range& param_min_max,
    const bool get_min_root, const bool try_limits_if_needed,
    PolynomialRoots& roots, bool invalid_solution, const MotionState& s) {
  double exec_time = 0.0;

  if (!invalid_solution) {
    // sort roots
    std::sort(roots.begin(), roots.end());
    // choose travel direction based on get_min_root flag
    int start_ind = 0;
    int ind_increment = 1;
    if (!get_min_root) {
      start_ind = roots.size() - 1;
      ind_increment = -1;
    }

    invalid_solution = true;
    int ind = start_ind;

    for (int i = 0; i < roots.size(); i++) {
      // check validity of solution
      if (IsInRange(param_min_max, roots[ind])) {
        // check if the min exec time is exceeded
        double tt = (texe_function)(s, roots[ind]);
        if (tt > (min_execution_time + kStep1TimeEpsilon)) {
          invalid_solution = false;
          exec_time = tt;
          break;
        }
      }
      ind += ind_increment;
    }
    // if no valid solution was found,
    // check whether any of the solutions is closest to one of the limits and
    // pick the closest one.
    if (invalid_solution) {
      ind = start_ind;
      double best_time = 0.0;
      double min_diff = GetRangeSize(param_min_max);
      for (int i = 0; i < roots.size(); i++) {
        // check validity of solution
        if ((roots[ind] >= param_min_max.min - kParamLimitOffset) &&
            (roots[ind] <= param_min_max.max + kParamLimitOffset)) {
          // check if the min exec time is exceeded
          double tt = (*texe_function)(s, roots[ind]);
          if (tt > (min_execution_time + kStep1TimeEpsilon)) {
            invalid_solution = false;
            double diff = std::min(fabs(roots[ind] - param_min_max.min),
                                   fabs(param_min_max.max - roots[ind]));
            if (diff < min_diff) {
              best_time = tt;
            }
          }
        }
        ind += ind_increment;
      }
      if (!invalid_solution) {
        exec_time = best_time;
      }
    }
  }
  // if no valid solution has been found yet,
  // check whether either of the limits works.
  if (invalid_solution && try_limits_if_needed) {
    double roots_new[2] = {param_min_max.min, param_min_max.max};
    double best_time = std::numeric_limits<double>::infinity();
    for (double i : roots_new) {
      // check if the min exec time is exceeded
      double tt = (*texe_function)(s, i);
      if (tt > (min_execution_time + kStep1TimeEpsilon)) {
        invalid_solution = false;
        if (tt < best_time) {
          best_time = tt;
        }
      }
    }
    if (!invalid_solution) {
      exec_time = best_time;
    }
  }

  if (invalid_solution) {
    return SetFailure(s);
  }

  return IncrementT(exec_time, s) | SetProfile(profile);
}

double GetExecutionTimeForPosTrapNegTrap(const MotionState& s,
                                         const double dt2) {
  return (Power2(GetAMax(s)) * (GetJMax(s) - GetJMin(s)) +
          Power2(GetAMin(s)) * (GetJMax(s) - GetJMin(s)) -
          2 * GetAMin(s) * (GetA(s) - dt2 * GetJMax(s)) * GetJMin(s) -
          2 * GetAMax(s) *
              (GetAMin(s) * (GetJMax(s) - GetJMin(s)) +
               dt2 * GetJMax(s) * GetJMin(s)) +
          GetJMin(s) * (Power2(GetA(s)) - 2 * GetJMax(s) * GetV(s) +
                        2 * GetJMax(s) * GetVTrgt(s))) /
         (2. * GetAMin(s) * GetJMax(s) * GetJMin(s));
}

double GetExecutionTimeForPosTrapNegTri(const MotionState& s,
                                        const double apeak) {
  return (0.5 * (2. * GetAMax(s) * apeak * GetJMax(s) -
                 2. * GetAMax(s) * (GetA(s) + apeak) * GetJMin(s) +
                 Power2(GetAMax(s)) * (-GetJMax(s) + GetJMin(s)) +
                 Power2(apeak) * (-GetJMax(s) + GetJMin(s)) +
                 GetJMin(s) * (Power2(GetA(s)) - 2. * GetJMax(s) * GetV(s) +
                               2. * GetJMax(s) * GetVTrgt(s)))) /
         (GetAMax(s) * GetJMax(s) * GetJMin(s));
}

double GetExecutionTimeForPosTriNegTrap(const MotionState& s,
                                        const double apeak) {
  return (Power2(GetAMin(s)) * (GetJMax(s) - GetJMin(s)) +
          Power2(apeak) * (GetJMax(s) - GetJMin(s)) -
          2 * GetAMin(s) *
              (apeak * (GetJMax(s) - GetJMin(s)) + GetA(s) * GetJMin(s)) +
          GetJMin(s) * (Power2(GetA(s)) - 2 * GetJMax(s) * GetV(s) +
                        2 * GetJMax(s) * GetVTrgt(s))) /
         (2. * GetAMin(s) * GetJMax(s) * GetJMin(s));
}

double GetExecutionTimeForNegLinPosTrap(const MotionState& s,
                                        const double aend) {
  return (Power2(GetA(s)) * GetJMax(s) - 2 * GetA(s) * GetAMax(s) * GetJMax(s) -
          Power2(GetAMax(s)) * GetJMax(s) +
          2 * aend * GetAMax(s) * (GetJMax(s) - GetJMin(s)) +
          Power2(GetAMax(s)) * GetJMin(s) +
          Power2(aend) * (-GetJMax(s) + GetJMin(s)) -
          2 * GetJMax(s) * GetJMin(s) * GetV(s) +
          2 * GetJMax(s) * GetJMin(s) * GetVTrgt(s)) /
         (2. * GetAMax(s) * GetJMax(s) * GetJMin(s));
}

Range GetPosTrapHoldTimeLimitsForPosTrapNegTrap(const MotionState& s) {
  double time_min, time_max;
  // if a->amax->0 and then a->amin->0, is v<=vtrgt?
  if (AUpToAMaxDownToZero(s) | ADownToAMinUpToZero | CmpVLteVTrgt) {
    // Profile: a -> amax -> hold -> 0 -> amin -> 0, so that v = vtrgt
    time_min = ((GetVTrgt(s) - GetV(s)) / GetAMax(s)) +
               (Power2(GetAMin(s)) * (GetJMax(s) - GetJMin(s)) /
                (2.0 * GetAMax(s) * GetJMax(s) * (-GetJMin(s)))) -
               (GetAMax(s) / (-2.0 * GetJMin(s))) -
               ((Power2(GetAMax(s)) - Power2(GetA(s))) /
                (2.0 * GetAMax(s) * GetJMax(s)));
  } else {
    time_min = 0.0;
  }
  time_min = std::max(time_min, 0.0);

  // Profile: a -> amax -> hold -> 0, so that v = vmax
  time_max = VUpToVMaxPosTrap(s) | GetT;

  time_max = std::max(time_max, 0.0);
  return {time_min, time_max};
}

// Limit Function for PosTrapNegTri Parameter: Peak Acceleration for NegTri
Range GetNegTriPeakAccelerationLimitsForPosTrapNegTri(const MotionState& s) {
  auto s_after_a_amax_zero = AUpToAMaxDownToZero(s);
  double pa_max = 0.0;
  if (CmpVGtVTrgt(s_after_a_amax_zero)) {
    // Profile: a -> amax -> 0 -> anegpeak1 -> 0, so that v = vtrgt
    pa_max = CalcPeakAccelVPrimeNegTri(GetVTrgt, s_after_a_amax_zero);
  }

  // Profile: v -> vmax (PosTrap), and then a -> anegpeak2 -> 0, so that v =
  // vtrgt
  auto s_after_vmax_postrap = SetV(GetVMax(s), s) | SetA(0);
  double pa_min = CalcPeakAccelVPrimeNegTri(GetVTrgt, s_after_vmax_postrap);

  return ClipToRange({pa_min, pa_max}, GetAMin(s), 0);
}

// Limit Function for PosTriNegTrap Parameter: Peak Acceleration for PosTri
Range GetPosTriPeakAccelerationLimitsForPosTriNegTrap(const MotionState& s) {
  double min = GetA(s);
  if (ADownToZero(s) | ADownToAMinUpToZero | CmpVGtVTrgt) {
    min = GetA(s);
  } else {
    // Profile: a -> apeak -> 0 -> amin -> 0, so that v = vtrgt
    min = GetSqrt(Power2(GetAMin(s)) + ((Power2(GetA(s)) / GetJMax(s)) +
                                        2.0 * (GetVTrgt(s) - GetV(s))) *
                                           (-GetJMin(s)) * GetJMax(s) /
                                           ((-GetJMin(s)) + GetJMax(s)));
  }

  // Ensure a <= peak accel min, max <= amax
  return ClipToRange({min, CalcPeakAccelVMaxPosTri(s)}, GetA(s), GetAMax(s));
}

Range GetNegLinEndAccelerationLimitsForNegLinPosTrap(const MotionState& s) {
  double min = 0.0;
  if (ADownToZero(s) | AUpToAMaxDownToZero | CmpVGteVTrgt) {
    // Profile: a -> aneglinend -> amax -> 0, so that v = vtrgt
    min =
        GetSqrt(Power2(GetAMax(s)) -
                ((2.0 * (GetVTrgt(s) - GetV(s)) * GetJMax(s) * (-GetJMin(s))) -
                 (Power2(GetA(s)) * GetJMax(s))) /
                    (GetJMax(s) + (-GetJMin(s))));
  }

  return ClipToRange({min, GetA(s)}, 0, GetA(s));
}

// Expand the range of parameters to valid solution from numerical solver if its
// within a certain absolute limit
inline Range ExpandIfWithinAbsLimitOffset(const Range& value_min_max) {
  // shrink the range of values
  if (fabs(GetRangeSize(value_min_max)) < kAbsoluteLimitOffset) {
    auto [value_min, value_max] = value_min_max;
    double offset_value =
        kRelativeLimitOffset * (value_max - value_min) + kAbsoluteLimitOffset;
    return {value_min - offset_value, value_max + offset_value};
  }

  return value_min_max;
}

}  // namespace

// Step1 Red Profile (PosTrapNegTrap)
// v->vtrgt and p->ptrgt
MotionState CalcExecTimePosTrapNegTrapFunc(const bool get_min_root,
                                           const bool try_limits_if_needed,
                                           const double min_execution_time,
                                           const MotionState& s) {
  // limits for param
  // param is hold time for postrap
  auto param_min_max = ExpandIfWithinAbsLimitOffset(
      GetPosTrapHoldTimeLimitsForPosTrapNegTrap(s));

  double pi = GetP(s);
  double vi = GetV(s);
  double ai = GetA(s);

  double amax = GetAMax(s);
  double jmax = GetJMax(s);

  double amin = GetAMin(s);
  double jmin = GetJMin(s);

  double vtrgt = GetVTrgt(s);
  double ptrgt = GetPTrgt(s);

  // postrap hold time
  std::array<double, 3> coefficients;
  coefficients[2] = amax * (0.5 - (0.5 * amax) / amin);
  coefficients[1] = (Power3(amax) * (0.5 * jmax - 0.5 * jmin) +
                     Power2(amax) * amin * (-jmax + 0.5 * jmin) +
                     amin * jmin * (-0.5 * Power2(ai) + jmax * vi) +
                     amax * (0.5 * Power2(amin) * jmax +
                             0.5 * Power2(ai) * jmin - jmax * jmin * vi)) /
                    (amin * jmax * jmin);
  coefficients[0] =
      (Power4(amax) *
           (-0.125 * Power2(jmax) + 0.25 * jmax * jmin - 0.125 * Power2(jmin)) +
       Power4(amin) *
           (kOneTwentyFourth * Power2(jmax) - kOneTwentyFourth * Power2(jmin)) +
       Power3(amax) * amin *
           (kOneThird * Power2(jmax) - 0.5 * jmax * jmin +
            kOneSixth * Power2(jmin)) +
       amin * Power2(jmin) *
           (kOneThird * Power3(ai) + Power2(jmax) * (pi - ptrgt) -
            ai * jmax * vi) +
       amax * amin * jmin *
           (Power2(ai) * (0.5 * jmax - 0.5 * jmin) +
            jmax * (-jmax + jmin) * vi) +
       Power2(amax) * (Power2(amin) * jmax * (-0.25 * jmax + 0.25 * jmin) +
                       jmin * (Power2(ai) * (-0.25 * jmax + 0.25 * jmin) +
                               jmax * (0.5 * jmax - 0.5 * jmin) * vi)) +
       Power2(amin) * jmax * jmin *
           (-0.25 * Power2(ai) + 0.5 * jmax * vi - 0.5 * jmin * vtrgt) +
       Power2(jmin) *
           (-0.125 * Power4(ai) + 0.5 * Power2(ai) * jmax * vi +
            Power2(jmax) * (-0.5 * Power2(vi) + 0.5 * Power2(vtrgt)))) /
      (amin * Power2(jmax) * Power2(jmin));

  // find roots
  auto roots = CalculatePolynomialRoots(coefficients);

  // find the execution time
  bool invalid_solution = (roots.size() != 2);
  return FindExecutionTime(Profile::kPosTrapNegTrap,
                           &GetExecutionTimeForPosTrapNegTrap,
                           min_execution_time, param_min_max, get_min_root,
                           try_limits_if_needed, roots, invalid_solution, s);
}

// Step1 Red Profile (PosTrapNegTri)
// v->vtrgt and p->ptrgt
MotionState CalcExecTimePosTrapNegTriFunc(const bool get_min_root,
                                          const bool try_limits_if_needed,
                                          const double min_execution_time,
                                          const MotionState& s) {
  // limits for param
  // param is (-ve) peak acceleration of negtri
  auto param_min_max = ExpandIfWithinAbsLimitOffset(
      GetNegTriPeakAccelerationLimitsForPosTrapNegTri(s));

  double pi = GetP(s);
  double vi = GetV(s);
  double ai = GetA(s);

  double amax = GetAMax(s);
  double jmax = GetJMax(s);

  double jmin = GetJMin(s);

  double vtrgt = GetVTrgt(s);
  double ptrgt = GetPTrgt(s);

  // negtri peak acceleration
  // quartic equation coefficients
  // a1*x^4 + b1*x^3 + c1*x^2 + d1*x + e1 = 0
  std::array<double, 5> coefficients;
  coefficients[4] =
      (0.125 * Power2(jmax) - 0.25 * jmax * jmin + 0.125 * Power2(jmin)) /
      (amax * Power2(jmax) * Power2(jmin));
  coefficients[3] = -kOneSixth / Power2(jmax) - kOneThird / Power2(jmin) +
                    0.5 / (jmax * jmin);
  coefficients[2] = (Power2(amax) * (0.25 * jmax - 0.25 * jmin) +
                     (-0.5 * jmax + 0.5 * jmin) * jmin * vtrgt) /
                    (amax * jmax * Power2(jmin));
  coefficients[1] = (-vtrgt) / jmax + (vtrgt) / jmin;
  coefficients[0] =
      (Power4(amax) * (-kOneTwentyFourth * Power2(jmax) +
                       kOneTwentyFourth * Power2(jmin)) +
       amax * Power2(jmin) *
           (kOneThird * Power3(ai) + Power2(jmax) * (pi - ptrgt) -
            ai * jmax * vi) +
       Power2(amax) * jmin *
           (-0.25 * Power2(ai) * jmin +
            jmax * (0.5 * jmin * vi - 0.5 * jmax * vtrgt)) +
       Power2(jmin) *
           (-0.125 * Power4(ai) + 0.5 * Power2(ai) * jmax * vi +
            Power2(jmax) * (-0.5 * Power2(vi) + 0.5 * Power2(vtrgt)))) /
      (amax * Power2(jmax) * Power2(jmin));

  // find roots
  PolynomialRoots roots = CalculatePolynomialRoots(coefficients);

  // find the execution time
  bool invalid_solution = false;
  if (roots.size() != 2 && roots.size() != 4) {
    invalid_solution = true;
  }

  return FindExecutionTime(Profile::kPosTrapNegTri,
                           &GetExecutionTimeForPosTrapNegTri,
                           min_execution_time, param_min_max, get_min_root,
                           try_limits_if_needed, roots, invalid_solution, s);
}

// Step1 Red Profile (PosTriNegTrap)
// v->vtrgt and p->ptrgt
MotionState CalcExecTimePosTriNegTrapFunc(const bool get_min_root,
                                          const bool try_limits_if_needed,
                                          const double min_execution_time,
                                          const MotionState& s) {
  // limits for param
  // param is peak acceleration of postri
  auto param_min_max = ExpandIfWithinAbsLimitOffset(
      GetPosTriPeakAccelerationLimitsForPosTriNegTrap(s));

  double pi = GetP(s);
  double vi = GetV(s);
  double ai = GetA(s);

  double jmax = GetJMax(s);

  double amin = GetAMin(s);
  double jmin = GetJMin(s);

  double vtrgt = GetVTrgt(s);
  double ptrgt = GetPTrgt(s);

  // postri peak acceleration
  // quartic equation coefficients
  // a1*x^4 + b1*x^3 + c1*x^2 + d1*x + e1 = 0
  std::array<double, 5> coefficients;
  coefficients[4] =
      (-0.125 * Power2(jmax) + 0.25 * jmax * jmin - 0.125 * Power2(jmin)) /
      (amin * Power2(jmax) * Power2(jmin));
  coefficients[3] =
      kOneSixth / Power2(jmax) + kOneThird / Power2(jmin) - 0.5 / (jmax * jmin);
  coefficients[2] = (Power2(amin) * jmax * (-0.25 * jmax + 0.25 * jmin) +
                     jmin * (Power2(ai) * (-0.25 * jmax + 0.25 * jmin) +
                             jmax * (0.5 * jmax - 0.5 * jmin) * vi)) /
                    (amin * Power2(jmax) * Power2(jmin));
  coefficients[1] =
      (Power2(ai) * (0.5 * jmax - 0.5 * jmin) + jmax * (-jmax + jmin) * vi) /
      (Power2(jmax) * jmin);
  coefficients[0] =
      (Power4(amin) *
           (kOneTwentyFourth * Power2(jmax) - kOneTwentyFourth * Power2(jmin)) +
       amin * Power2(jmin) *
           (kOneThird * Power3(ai) + Power2(jmax) * (pi - ptrgt) -
            ai * jmax * vi) +
       Power2(amin) * jmax * jmin *
           (-0.25 * Power2(ai) + 0.5 * jmax * vi - 0.5 * jmin * vtrgt) +
       Power2(jmin) *
           (-0.125 * Power4(ai) + 0.5 * Power2(ai) * jmax * vi +
            Power2(jmax) * (-0.5 * Power2(vi) + 0.5 * Power2(vtrgt)))) /
      (amin * Power2(jmax) * Power2(jmin));

  // find roots
  PolynomialRoots roots = CalculatePolynomialRoots(coefficients);

  // find the execution time
  bool invalid_solution = false;
  if (roots.size() != 2 && roots.size() != 4) {
    invalid_solution = true;
  }
  return FindExecutionTime(Profile::kPosTriNegTrap,
                           &GetExecutionTimeForPosTriNegTrap,
                           min_execution_time, param_min_max, get_min_root,
                           try_limits_if_needed, roots, invalid_solution, s);
}

// Step1 Red Profile (NegLinPosTrap)
// v->vtrgt and p->ptrgt
MotionState CalcExecTimeNegLinPosTrapFunc(const bool get_min_root,
                                          const bool try_limits_if_needed,
                                          const double min_execution_time,
                                          const MotionState& s) {
  // limits for param
  // param is end acceleration of neglin
  auto param_min_max = ExpandIfWithinAbsLimitOffset(
      GetNegLinEndAccelerationLimitsForNegLinPosTrap(s));

  double pi = GetP(s);
  double vi = GetV(s);
  double ai = GetA(s);

  double amax = GetAMax(s);

  double jmin = GetJMin(s);
  double jmax = GetJMax(s);

  double vtrgt = GetVTrgt(s);
  double ptrgt = GetPTrgt(s);

  // neglin end acceleration
  // quartic equation coefficients
  // a1*x^4 + b1*x^3 + c1*x^2 + d1*x + e1 = 0
  std::array<double, 5> coefficients;
  coefficients[4] =
      (-0.125 * Power2(jmax) + 0.25 * jmax * jmin - 0.125 * Power2(jmin)) /
      (amax * Power2(jmax) * Power2(jmin));
  coefficients[3] =
      kOneThird / Power2(jmax) + kOneSixth / Power2(jmin) - 0.5 / (jmax * jmin);
  coefficients[2] = (Power2(ai) * jmax * (0.25 * jmax - 0.25 * jmin) +
                     jmin * (Power2(amax) * (0.25 * jmax - 0.25 * jmin) +
                             jmax * (-0.5 * jmax + 0.5 * jmin) * vi)) /
                    (amax * Power2(jmax) * Power2(jmin));
  coefficients[1] =
      (Power2(ai) * (-0.5 * jmax + 0.5 * jmin) + (jmax - jmin) * jmin * vi) /
      (jmax * Power2(jmin));
  coefficients[0] =
      Power3(amax) *
          (kOneTwentyFourth / Power2(jmax) - kOneTwentyFourth / Power2(jmin)) +
      pi - ptrgt + (kOneThird * Power3(ai) - ai * jmin * vi) / Power2(jmin) +
      (amax * (-0.25 * Power2(ai) + 0.5 * jmin * vi - 0.5 * jmax * vtrgt)) /
          (jmax * jmin) +
      ((-0.125 * Power4(ai)) / Power2(jmin) + (0.5 * Power2(ai) * vi) / jmin -
       0.5 * Power2(vi) + 0.5 * Power2(vtrgt)) /
          amax;

  // find roots
  PolynomialRoots roots = CalculatePolynomialRoots(coefficients);

  // find the execution time
  bool invalid_solution = false;
  if (roots.size() != 2 && roots.size() != 4) {
    invalid_solution = true;
  }

  return FindExecutionTime(Profile::kNegLinPosTrap,
                           &GetExecutionTimeForNegLinPosTrap,
                           min_execution_time, param_min_max, get_min_root,
                           try_limits_if_needed, roots, invalid_solution, s);
}
Range GetPosTriPeakAccelerationLimitsForPosTriNegTri(const MotionState& s) {
  Range result;

  if (ADownToZero(s) | CmpVGtVTrgt) {
    result.min = GetA(s);
  } else {
    // Profile: a -> apeak -> 0, so that v = vtrgt
    result.min = CalcPeakAccelVTrgtPosTri(s);
  }

  // Profile: a -> apeak -> 0 -> amin -> 0, so that v = vtrgt
  double postri_peak_acceleration_maxduetoamax =
      GetSqrt(Power2(GetAMin(s)) +
              (2.0 * (GetVTrgt(s) - GetV(s)) + (Power2(GetA(s)) / GetJMax(s))) *
                  (-GetJMin(s)) * GetJMax(s) / ((-GetJMin(s)) + GetJMax(s)));

  postri_peak_acceleration_maxduetoamax =
      std::min(postri_peak_acceleration_maxduetoamax, GetAMax(s));

  // Profile: a -> apeak -> 0, so that v = vmax
  double postri_peak_acceleration_maxduetovmax = CalcPeakAccelVMaxPosTri(s);

  postri_peak_acceleration_maxduetovmax =
      std::min(postri_peak_acceleration_maxduetovmax, GetAMax(s));

  result.max = std::min(postri_peak_acceleration_maxduetovmax,
                        postri_peak_acceleration_maxduetoamax);

  return ClipToRange(result, GetA(s), GetAMax(s));
}

Range GetNegLinEndAccelerationLimitsForNegLinPosTri(const MotionState& s) {
  double max;
  if (AUpToAMaxDownToZero(s) | CmpVGtVTrgt) {
    max = GetA(s);
  } else {
    // Profile: a -> aneglinend -> amax -> 0, so that v = vtrgt
    max =
        GetSqrt(Power2(GetAMax(s)) -
                ((2.0 * (GetVTrgt(s) - GetV(s)) * GetJMax(s) * (-GetJMin(s))) -
                 (Power2(GetA(s)) * GetJMax(s))) /
                    (GetJMax(s) + (-GetJMin(s))));
  }

  return ClipToRange({0, max}, 0, GetA(s));
}

double GetPositionErrorForPosTriNegTriFunc(
    const MotionState& s, const double peak_acceleration_postri,
    bool* invalid_solution, double* execution_time) {
  auto s_after_postri =
      SetT(0, s) | AUpToAPrimePosTri(peak_acceleration_postri);

  // check whether the solution is valid
  if (invalid_solution != nullptr) {
    if ((GetV(s_after_postri) - GetVMax(s)) > kAbsVelocityErrorTolerance) {
      (*invalid_solution) = true;
    } else {
      (*invalid_solution) = false;
    }
  }

  auto s_after_postri_capped = s_after_postri | CapVToVMax;

  // NegTri

  // check whether the solution is valid
  if (invalid_solution != nullptr) {
    if ((GetV(s_after_postri_capped) - GetVTrgt(s)) <
        -kAbsVelocityErrorTolerance) {
      (*invalid_solution) = true;
    } else {
      (*invalid_solution) = false;
    }
  }

  auto s_after_neg_tri = s_after_postri_capped | VDownToVTrgtNegTri;

  // return execution time
  if (execution_time != nullptr) {
    (*execution_time) = GetT(s_after_neg_tri);
  }
  // return position error
  return (GetP(s_after_neg_tri) - GetPTrgt(s));
}

double GetDerivativePositionErrorForPosTriNegTriFunc(
    const MotionState& s, const double peak_acceleration_postri,
    bool* invalid_solution, double* execution_time) {
  // PosTri

  // change in velocity during acceleration and deceleration phases of PosTri
  double delta_velocity_acceleration_phase_postri =
      (Power2(peak_acceleration_postri) - Power2(GetA(s))) * GetHalfInvJMax(s);
  double delta_velocity_deceleration_phase_postri =
      (Power2(peak_acceleration_postri) * GetNegHalfInvJMin(s));

  // partial derivative of the above change in velocities of PosTri w.r.t.
  // peak_acceleration_postri
  double
      d_delta_velocity_acceleration_phase_postri__d_peak_acceleration_postri =
          peak_acceleration_postri * GetInvJMax(s);
  double
      d_delta_velocity_deceleration_phase_postri__d_peak_acceleration_postri =
          peak_acceleration_postri * GetNegInvJMin(s);

  // time durations of acceleration and deceleration phases of PosTri
  double acceleration_time_postri =
      (peak_acceleration_postri - GetA(s)) * GetInvJMax(s);
  double deceleration_time_postri = peak_acceleration_postri * GetNegInvJMin(s);

  // time duration of PosTri
  double delta_time_postri =
      acceleration_time_postri + deceleration_time_postri;

  // resulting position, velocity and acceleration after acceleration phase of
  // PosTri and their partial derivatives w.r.t. peak_acceleration_postri (which
  // is fixed at inv_jmax for accel and -inv_jmin for decel)
  double d_resulting_position__d_peak_acceleration_postri =
      (GetV(s) + GetA(s) * acceleration_time_postri +
       GetHalfJMax(s) * Power2(acceleration_time_postri)) *
      GetInvJMax(s);

  double resulting_velocity =
      GetV(s) + delta_velocity_acceleration_phase_postri;
  double d_resulting_velocity__d_peak_acceleration_postri =
      d_delta_velocity_acceleration_phase_postri__d_peak_acceleration_postri;
  d_resulting_position__d_peak_acceleration_postri +=
      (resulting_velocity +
       (peak_acceleration_postri * deceleration_time_postri) +
       (GetHalfJMin(s) * Power2(deceleration_time_postri))) *
          GetNegInvJMin(s) +
      ((0.5 * deceleration_time_postri) +
       d_resulting_velocity__d_peak_acceleration_postri) *
          deceleration_time_postri;

  resulting_velocity += delta_velocity_deceleration_phase_postri;
  d_resulting_velocity__d_peak_acceleration_postri +=
      d_delta_velocity_deceleration_phase_postri__d_peak_acceleration_postri;

  // check whether the solution is valid
  if (invalid_solution != nullptr) {
    if ((resulting_velocity - GetVMax(s)) > kAbsVelocityErrorTolerance) {
      (*invalid_solution) = true;
    } else {
      (*invalid_solution) = false;
    }
  }
  if (resulting_velocity > GetVMax(s)) {
    resulting_velocity = GetVMax(s);
  }

  // NegTri

  // check whether the solution is valid
  if (invalid_solution != nullptr) {
    if ((resulting_velocity - GetVTrgt(s)) < -kAbsVelocityErrorTolerance) {
      (*invalid_solution) = true;
    } else {
      (*invalid_solution) = false;
    }
  }

  // peak_acceleration_negtri needed for NegTri such that v=vtrgt after NegTri
  // and its partial derivative w.r.t. peak_acceleration_postri
  double peak_acceleration_negtri = 0.0;
  double d_peak_acceleration_negtri__d_peak_acceleration_postri = 0.0;

  if (resulting_velocity > GetVTrgt(s)) {
    peak_acceleration_negtri =
        -1.0 *
        GetSqrt(2.0 * (resulting_velocity - GetVTrgt(s)) * GetJerkConstant(s));
    d_peak_acceleration_negtri__d_peak_acceleration_postri =
        -1.0 *
        GetSqrt(GetJerkConstant(s) /
                (2.0 * (resulting_velocity - GetVTrgt(s)))) *
        d_resulting_velocity__d_peak_acceleration_postri;
  }

  // change in velocity during deceleration and acceleration phases of NegTri
  double delta_velocity_deceleration_phase_negtri =
      -(Power2(peak_acceleration_negtri) * GetNegHalfInvJMin(s));

  // partial derivative of the above change in velocities of NegTri w.r.t.
  // peak_acceleration_postri
  double
      d_delta_velocity_deceleration_phase_negtri__d_peak_acceleration_postri =
          (peak_acceleration_negtri *
           d_peak_acceleration_negtri__d_peak_acceleration_postri) *
          GetInvJMin(s);

  // time durations of deceleration, hold and acceleration phases of NegTri
  double deceleration_time_negtri = peak_acceleration_negtri * GetInvJMin(s);
  double acceleration_time_negtri = (-peak_acceleration_negtri) * GetInvJMax(s);

  // time duration of NegTri
  double delta_time_negtri =
      deceleration_time_negtri + acceleration_time_negtri;

  // partial derivative of time durations of deceleration, hold and acceleration
  // phases of NegTri w.r.t. peak_acceleration_postri
  double d_deceleration_time_negtri__d_peak_acceleration_postri =
      d_peak_acceleration_negtri__d_peak_acceleration_postri * GetInvJMin(s);
  double d_acceleration_time_negtri__d_peak_acceleration_postri =
      -d_peak_acceleration_negtri__d_peak_acceleration_postri * GetInvJMax(s);

  // resulting position, velocity and acceleration after deceleration phase of
  // NegTri and their partial derivatives w.r.t. peak_acceleration_postri
  d_resulting_position__d_peak_acceleration_postri +=
      (resulting_velocity +
       (GetHalfJMin(s) * Power2(deceleration_time_negtri))) *
          d_deceleration_time_negtri__d_peak_acceleration_postri +
      (d_resulting_velocity__d_peak_acceleration_postri *
       deceleration_time_negtri);

  resulting_velocity += delta_velocity_deceleration_phase_negtri;
  d_resulting_velocity__d_peak_acceleration_postri +=
      d_delta_velocity_deceleration_phase_negtri__d_peak_acceleration_postri;

  // partial derivative of resulting position after acceleration phase of NegTri
  // w.r.t peak_acceleration_postri
  d_resulting_position__d_peak_acceleration_postri +=
      (resulting_velocity +
       (peak_acceleration_negtri * acceleration_time_negtri) +
       (GetHalfJMax(s) * Power2(acceleration_time_negtri))) *
          d_acceleration_time_negtri__d_peak_acceleration_postri +
      acceleration_time_negtri *
          (d_resulting_velocity__d_peak_acceleration_postri +
           (0.5 * acceleration_time_negtri *
            d_peak_acceleration_negtri__d_peak_acceleration_postri));

  // return execution time
  if (execution_time != nullptr) {
    (*execution_time) = delta_time_postri + delta_time_negtri;
  }

  // return derivative of position error
  return d_resulting_position__d_peak_acceleration_postri;
}

double GetPositionErrorForNegLinPosTriFuncFunc(
    const MotionState& s, const double end_acceleration_neglin,
    bool* invalid_solution, double* execution_time) {
  auto s1 = ADownToAPrime(end_acceleration_neglin, s);

  // check whether the solution is valid
  if (invalid_solution != nullptr) {
    if ((GetVTrgt(s) - GetV(s1)) < -kAbsVelocityErrorTolerance) {
      (*invalid_solution) = true;
    } else {
      (*invalid_solution) = false;
    }
  }

  auto s2 = VUpToVTrgtPosTri(s1);

  // return execution time
  if (execution_time != nullptr) {
    (*execution_time) = GetT(s2) - GetT(s);
  }

  // return position error
  return (GetP(s2) - GetPTrgt(s));
}

double GetDerivativePositionErrorForNegLinPosTriFunc(
    const MotionState& s, const double end_acceleration_neglin,
    bool* invalid_solution, double* execution_time) {
  // NegLin

  // time duration of Neglin
  double delta_time_neglin =
      (GetA(s) - end_acceleration_neglin) * GetNegInvJMin(s);

  // resulting position, velocity and acceleration after NegLin and
  // their partial derivatives w.r.t. end_acceleration_neglin (constant
  // inv_jmin)
  double d_resulting_position__d_end_acceleration_neglin =
      (GetV(s) + GetA(s) * delta_time_neglin +
       GetHalfJMin(s) * Power2(delta_time_neglin)) *
      GetInvJMin(s);

  double resulting_velocity = GetV(s) + GetA(s) * delta_time_neglin +
                              GetHalfJMin(s) * Power2(delta_time_neglin);

  double d_resulting_velocity__d_end_acceleration_neglin =
      end_acceleration_neglin / GetJMin(s);
  // PosTri

  // check whether the solution is valid
  if (invalid_solution != nullptr) {
    if ((GetVTrgt(s) - resulting_velocity) < -kAbsVelocityErrorTolerance) {
      (*invalid_solution) = true;
    } else {
      (*invalid_solution) = false;
    }
  }

  // peak_acceleration_postri needed for PosTri such that v=vtrgt after PosTri
  // and partial derivative of peak_acceleration_postri w.r.t.
  // peak_acceleration_postri
  double peak_acceleration_postri = 0.0;
  double d_peak_acceleration_postri__d_end_acceleration_neglin = 0.0;
  if (resulting_velocity < GetVTrgt(s)) {
    peak_acceleration_postri =
        GetSqrt((2.0 * (GetVTrgt(s) - resulting_velocity) +
                 (Power2(end_acceleration_neglin) * GetInvJMax(s))) *
                GetJerkConstant(s));

    d_peak_acceleration_postri__d_end_acceleration_neglin =
        end_acceleration_neglin / peak_acceleration_postri;
  }

  // change in velocity during acceleration and deceleration phases of PosTri
  double delta_velocity_acceleration_phase_postri =
      (Power2(peak_acceleration_postri) - Power2(end_acceleration_neglin)) *
      GetHalfInvJMax(s);

  // partial derivative of the above change in velocity w.r.t.
  // end_acceleration_neglin
  double d_delta_velocity_acceleration_phase_postri__d_end_acceleration_neglin =
      ((peak_acceleration_postri *
        d_peak_acceleration_postri__d_end_acceleration_neglin) -
       end_acceleration_neglin) *
      GetInvJMax(s);

  // time durations of acceleration and deceleration phases of PosTri
  double acceleration_time_postri =
      (peak_acceleration_postri - end_acceleration_neglin) * GetInvJMax(s);
  double deceleration_time_postri = peak_acceleration_postri * GetNegInvJMin(s);

  // time duration of PosTri
  double delta_time_postri =
      acceleration_time_postri + deceleration_time_postri;

  // partial derivative of the above time durations w.r.t.
  // end_acceleration_neglin
  double d_acceleration_time_postri__d_end_acceleration_neglin =
      (d_peak_acceleration_postri__d_end_acceleration_neglin - 1.0) *
      GetInvJMax(s);
  double d_deceleration_time_postri__d_end_acceleration_neglin =
      d_peak_acceleration_postri__d_end_acceleration_neglin * GetNegInvJMin(s);

  // resulting position, velocity and acceleration after acceleration phase of
  // PosTri and their partial derivatives w.r.t. end_acceleration_neglin
  d_resulting_position__d_end_acceleration_neglin +=
      (resulting_velocity +
       (end_acceleration_neglin * acceleration_time_postri) +
       GetHalfJMax(s) * Power2(acceleration_time_postri)) *
          d_acceleration_time_postri__d_end_acceleration_neglin +
      ((0.5 * acceleration_time_postri) +
       d_resulting_velocity__d_end_acceleration_neglin) *
          acceleration_time_postri;

  resulting_velocity += delta_velocity_acceleration_phase_postri;
  d_resulting_velocity__d_end_acceleration_neglin +=
      d_delta_velocity_acceleration_phase_postri__d_end_acceleration_neglin;

  // partial derivative of resulting position, velocity and acceleration after
  // deceleration phase of PosTri w.r.t. end_acceleration_neglin
  d_resulting_position__d_end_acceleration_neglin +=
      (resulting_velocity +
       (peak_acceleration_postri * deceleration_time_postri) +
       (GetHalfJMin(s) * Power2(deceleration_time_postri))) *
          d_deceleration_time_postri__d_end_acceleration_neglin +
      deceleration_time_postri *
          (d_resulting_velocity__d_end_acceleration_neglin +
           (0.5 * d_peak_acceleration_postri__d_end_acceleration_neglin *
            deceleration_time_postri));

  // return execution time
  if (execution_time != nullptr) {
    (*execution_time) = delta_time_neglin + delta_time_postri;
  }

  // return derivative of position error
  return d_resulting_position__d_end_acceleration_neglin;
}

double GetSecondDerivativePositionErrorForNegLinPosTriFunc(
    const MotionState& s, const double end_acceleration_neglin,
    bool* invalid_solution, double* execution_time) {
  // NegLin

  // time duration of Neglin
  double delta_time_neglin =
      (GetA(s) - end_acceleration_neglin) * GetNegInvJMin(s);

  // resulting position, velocity and acceleration after NegLin and
  // their partial derivatives w.r.t. end_acceleration_neglin (i.e. inv_jmin)
  double dd_resulting_position__dd_end_acceleration_neglin =
      (GetA(s) + GetJMin(s) * delta_time_neglin) * Power2(GetInvJMin(s));

  double resulting_velocity = GetV(s) + GetA(s) * delta_time_neglin +
                              GetHalfJMin(s) * Power2(delta_time_neglin);
  double d_resulting_velocity__d_end_acceleration_neglin =
      end_acceleration_neglin * GetInvJMin(s);
  double dd_resulting_velocity__dd_end_acceleration_neglin = GetInvJMin(s);

  // PosTri

  // check whether the solution is valid
  if (invalid_solution != nullptr) {
    if ((GetVTrgt(s) - resulting_velocity) < -kAbsVelocityErrorTolerance) {
      (*invalid_solution) = true;
    } else {
      (*invalid_solution) = false;
    }
  }

  // peak_acceleration_postri needed for PosTri such that v=vtrgt after PosTri,
  // partial derivative of peak_acceleration_postri w.r.t.
  // peak_acceleration_postri and second partial derivative of
  // peak_acceleration_postri w.r.t. peak_acceleration_postri
  double peak_acceleration_postri = 0.0;
  double d_peak_acceleration_postri__d_end_acceleration_neglin = 0.0;
  double dd_peak_acceleration_postri__dd_end_acceleration_neglin = 0.0;
  if (resulting_velocity < GetVTrgt(s)) {
    peak_acceleration_postri =
        GetSqrt((2.0 * (GetVTrgt(s) - resulting_velocity) +
                 (Power2(end_acceleration_neglin) * GetInvJMax(s))) *
                GetJerkConstant(s));
    d_peak_acceleration_postri__d_end_acceleration_neglin =
        end_acceleration_neglin / peak_acceleration_postri;
    dd_peak_acceleration_postri__dd_end_acceleration_neglin =
        (1.0 / peak_acceleration_postri) -
        (Power2(end_acceleration_neglin) / Power3(peak_acceleration_postri));
  }

  // change in velocity during acceleration and deceleration phases of PosTri
  double delta_velocity_acceleration_phase_postri =
      (Power2(peak_acceleration_postri) - Power2(end_acceleration_neglin)) *
      GetHalfInvJMax(s);

  // partial derivative of the above change in velocity w.r.t.
  // end_acceleration_neglin
  double d_delta_velocity_acceleration_phase_postri__d_end_acceleration_neglin =
      ((peak_acceleration_postri *
        d_peak_acceleration_postri__d_end_acceleration_neglin) -
       end_acceleration_neglin) *
      GetInvJMax(s);

  // second partial derivative of the above change in velocity w.r.t.
  // end_acceleration_neglin
  double
      dd_delta_velocity_acceleration_phase_postri__dd_end_acceleration_neglin =
          ((peak_acceleration_postri *
            dd_peak_acceleration_postri__dd_end_acceleration_neglin) +
           Power2(d_peak_acceleration_postri__d_end_acceleration_neglin) -
           1.0) *
          GetInvJMax(s);

  // time durations of acceleration and deceleration phases of PosTri
  double acceleration_time_postri =
      (peak_acceleration_postri - end_acceleration_neglin) * GetInvJMax(s);
  double deceleration_time_postri =
      (peak_acceleration_postri * GetNegInvJMin(s));

  // time duration of PosTri
  double delta_time_postri =
      acceleration_time_postri + deceleration_time_postri;

  // partial derivative of the above time durations w.r.t.
  // end_acceleration_neglin
  double d_acceleration_time_postri__d_end_acceleration_neglin =
      (d_peak_acceleration_postri__d_end_acceleration_neglin - 1.0) *
      GetInvJMax(s);
  double d_deceleration_time_postri__d_end_acceleration_neglin =
      d_peak_acceleration_postri__d_end_acceleration_neglin * GetNegInvJMin(s);

  // second partial derivative of the above time durations w.r.t.
  // end_acceleration_neglin
  double dd_acceleration_time_postri__dd_end_acceleration_neglin =
      dd_peak_acceleration_postri__dd_end_acceleration_neglin * GetInvJMax(s);
  double dd_deceleration_time_postri__dd_end_acceleration_neglin =
      dd_peak_acceleration_postri__dd_end_acceleration_neglin *
      GetNegInvJMin(s);

  // resulting position, velocity and acceleration after acceleration phase of
  // PosTri and their partial derivatives w.r.t. end_acceleration_neglin

  dd_resulting_position__dd_end_acceleration_neglin +=
      dd_acceleration_time_postri__dd_end_acceleration_neglin *
          (resulting_velocity +
           (end_acceleration_neglin * acceleration_time_postri) +
           (GetHalfJMax(s) * Power2(acceleration_time_postri))) +
      d_acceleration_time_postri__d_end_acceleration_neglin *
          ((2.0 * (d_resulting_velocity__d_end_acceleration_neglin +
                   acceleration_time_postri)) +
           (d_acceleration_time_postri__d_end_acceleration_neglin *
            (end_acceleration_neglin +
             (GetJMax(s) * acceleration_time_postri)))) +
      dd_resulting_velocity__dd_end_acceleration_neglin *
          acceleration_time_postri;

  resulting_velocity += delta_velocity_acceleration_phase_postri;
  d_resulting_velocity__d_end_acceleration_neglin +=
      d_delta_velocity_acceleration_phase_postri__d_end_acceleration_neglin;
  dd_resulting_velocity__dd_end_acceleration_neglin +=
      dd_delta_velocity_acceleration_phase_postri__dd_end_acceleration_neglin;

  // second partial derivative of resulting position after deceleration phase
  // of PosTri w.r.t. end_acceleration_neglin
  dd_resulting_position__dd_end_acceleration_neglin +=
      resulting_velocity *
          dd_deceleration_time_postri__dd_end_acceleration_neglin +
      2.0 * d_resulting_velocity__d_end_acceleration_neglin *
          d_deceleration_time_postri__d_end_acceleration_neglin +
      dd_resulting_velocity__dd_end_acceleration_neglin *
          deceleration_time_postri +
      peak_acceleration_postri * deceleration_time_postri *
          dd_deceleration_time_postri__dd_end_acceleration_neglin +
      peak_acceleration_postri *
          Power2(d_deceleration_time_postri__d_end_acceleration_neglin) +
      2.0 * d_peak_acceleration_postri__d_end_acceleration_neglin *
          deceleration_time_postri *
          d_deceleration_time_postri__d_end_acceleration_neglin +
      0.5 * dd_peak_acceleration_postri__dd_end_acceleration_neglin *
          Power2(deceleration_time_postri) +
      GetJMin(s) * deceleration_time_postri *
          Power2(d_deceleration_time_postri__d_end_acceleration_neglin) +
      GetHalfJMin(s) * Power2(deceleration_time_postri) *
          dd_deceleration_time_postri__dd_end_acceleration_neglin;

  // return execution time
  if (execution_time != nullptr) {
    (*execution_time) = delta_time_neglin + delta_time_postri;
  }

  // return second derivative of position error
  return dd_resulting_position__dd_end_acceleration_neglin;
}

// Provides an additional error tolerance value in order to improve the
// numerical robustness
// Returns the absolute value of the additionally allowed position error.
double GetAnAdditionalNumericalErrorToleranceBasedOnTheExecutionTime(
    const double execution_time_value) {
  double return_value = 1.0, time_value = execution_time_value;

  if (time_value > kLargeExecutionTimeLevel1) {
    time_value *= kLargeExecutionTimeReductionRatio;

    for (unsigned int i = 0;
         i < kAdditionalPositionErrorToleranceRevisionIterations; i++) {
      if (time_value > kLargeExecutionTimeLevel2) {
        time_value *= kLargeExecutionTimeReductionRatio;
        return_value *= kAdditionalPositionErrorToleranceIncreaseRatio;
      } else {
        break;
      }
    }

    return return_value;
  } else {
    return 0.0;
  }
}

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic
