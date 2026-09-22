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

#include "intrinsic/icon/reflexxes/internal/position_calc_traj_funcs.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <utility>

#include "absl/functional/bind_front.h"
#include "intrinsic/icon/reflexxes/internal/constants.h"
#include "intrinsic/icon/reflexxes/internal/math.h"
#include "intrinsic/icon/reflexxes/internal/motion_state.h"
#include "intrinsic/icon/reflexxes/internal/motion_state_transforms.h"
#include "intrinsic/icon/reflexxes/internal/motion_state_with_polys.h"
#include "intrinsic/icon/reflexxes/internal/polynomial_solvers.h"
#include "intrinsic/icon/reflexxes/internal/synchronize_dofs_base.h"
#include "intrinsic/icon/reflexxes/profile.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {
namespace {

// The (integer) number of maximal loops for the
// Anderson-Bjoerck-King method within the blue profiles.
constexpr unsigned int kNumberOfLoopsABKBlue = 60;

// Return the difference between the current and target position, its error.
inline double GetPError(const SyncDOFsMotionState& s) {
  return GetP(s) - GetPTrgt(s);
}

inline Range ClipAndVerifyRange(const Range& range, const Range& max_limits) {
  Range clipped = ClipToRange(range, max_limits);
  clipped.max = std::max(clipped.min, clipped.max);
  return clipped;
}

// Limit Function for PosTrapZeroNegTrap
Range GetPosTrapHoldTimeLimitsForPosTrapZeroNegTrap(
    const double desired_execution_time, const SyncDOFsMotionState& state) {
  double min = 0., max = 0.;

  // The maximum hold time value can be determined by calculating the hold time
  // of PosTrapNegTrap such that t=tsync and v=vtrgt
  max = (GetVTrgt(state) - GetV(state) -
         0.5 * GetInvJMaxSubInvJMin(state) *
             (Power2(GetAMax(state)) - Power2(GetAMin(state))) +
         0.5 * Power2(GetA(state)) / GetJMax(state) -
         GetAMin(state) *
             (desired_execution_time -
              GetInvJMaxSubInvJMin(state) * (GetAMax(state) - GetAMin(state)) +
              GetA(state) / GetJMax(state))) /
        (GetAMax(state) - GetAMin(state));

  // The hold time required to achieve the desired execution time without the
  // vtrgt constraint using PosTrap and a(=0)->amin->0
  double hold_time_postrap_without_vtrgt =
      desired_execution_time -
      GetInvJMaxSubInvJMin(state) * (GetAMax(state) - GetAMin(state)) +
      GetA(state) / GetJMax(state);
  max = std::min(max, hold_time_postrap_without_vtrgt);

  // If vtrgt cannot be reached by a->amax->0->amin->0,
  // then the hold time of PosTrap is greater than zero
  double delta_velocity =
      0.5 * GetInvJMaxSubInvJMin(state) *
          (Power2(GetAMax(state)) - Power2(GetAMin(state))) -
      0.5 * Power2(GetA(state)) / GetJMax(state);

  if (delta_velocity <= (GetVTrgt(state) - GetV(state))) {
    min = (GetVTrgt(state) - GetV(state) - delta_velocity) / GetAMax(state);
  }

  return ClipAndVerifyRange({min, max}, {0., desired_execution_time});
}

Range GetNegTriPeakAccelerationLimitsForPosTrapZeroNegTri(
    const double desired_execution_time, const SyncDOFsMotionState& s) {
  double min = 0., max = 0.;

  // Profile: a -> amax -> 0
  double end_velocity_a2amax2zero =
      GetV(s) + (Power2(GetAMax(s)) - Power2(GetA(s))) / (2.0 * GetJMax(s)) +
      Power2(GetAMax(s)) / (-2.0 * GetJMin(s));

  if (end_velocity_a2amax2zero > GetVTrgt(s)) {
    // Profile: a -> amax -> 0 -> anegpeak1 -> 0, so that v = vtrgt
    max = -GetSqrt(2 * GetJMax(s) * (-GetJMin(s)) *
                   (end_velocity_a2amax2zero - GetVTrgt(s)) /
                   (GetJMax(s) + (-GetJMin(s))));
  }

  // The minimum (max in magnitude) peak acceleration value can be determined by
  // calculating the peak acceleration of PosTrapNegTri such that t=tsync and
  // v=vtrgt constant value
  std::array<double, 3> coefficients;
  coefficients[2] = -0.5 * GetInvJMaxSubInvJMin(s);
  coefficients[1] = GetAMax(s) * GetInvJMaxSubInvJMin(s);
  coefficients[0] = GetAMax(s) * (desired_execution_time -
                                  0.5 * GetInvJMaxSubInvJMin(s) * GetAMax(s) +
                                  GetA(s) / GetJMax(s)) -
                    (GetVTrgt(s) - GetV(s)) -
                    0.5 * Power2(GetA(s)) / GetJMax(s);

  PolynomialRoots roots = SafeCalculatePolynomialRoots(coefficients);
  if (roots[0] < 0.0) {
    min = roots[0];
  } else {
    if (roots[1] < 0.0) {
      min = roots[1];
    }
  }

  // The peak acceleration required to achieve the desired execution time
  // without the vtrgt constraint using a->amax->0 and NegTri
  double negtri_peak_acceleration_without_vtrgt =
      -(desired_execution_time - GetAMax(s) * GetInvJMaxSubInvJMin(s) +
        GetA(s) / GetJMax(s)) /
      GetInvJMaxSubInvJMin(s);
  min = std::max(min, negtri_peak_acceleration_without_vtrgt);

  return ClipAndVerifyRange({min, max}, {GetAMin(s), 0.0});
}

Range GetPosTriPeakAccelerationLimitsForPosTriZeroNegTrap(
    const double desired_execution_time, const SyncDOFsMotionState& s) {
  double min = 0., max = 0.;

  // a->0 and then a(=0)->amin->0
  double end_velocity_a2zero2amin2zero =
      GetV(s) + (Power2(GetA(s)) / (-2.0 * GetJMin(s))) -
      0.5 * Power2(GetAMin(s)) *
          ((GetJMax(s) + (-GetJMin(s))) / (GetJMax(s) * (-GetJMin(s))));

  if (end_velocity_a2zero2amin2zero > GetVTrgt(s)) {
    min = GetA(s);
  } else {
    // Profile: a -> apeak -> 0 -> amin -> 0, so that v = vtrgt
    min = GetSqrt(Power2(GetAMin(s)) + ((Power2(GetA(s)) / GetJMax(s)) +
                                        2.0 * (GetVTrgt(s) - GetV(s))) *
                                           (-GetJMin(s)) * GetJMax(s) /
                                           ((-GetJMin(s)) + GetJMax(s)));
  }

  // The maximum peak acceleration value can be determined by calculating
  // the peak acceleration of PosTriNegTrap such that t=tsync and v=vtrgt
  // constant value

  std::array<double, 3> coefficients;
  coefficients[2] = 0.5 * GetInvJMaxSubInvJMin(s);
  coefficients[1] = -GetAMin(s) * GetInvJMaxSubInvJMin(s);
  coefficients[0] = GetAMin(s) * (desired_execution_time +
                                  GetInvJMaxSubInvJMin(s) * GetAMin(s) +
                                  GetA(s) / GetJMax(s)) -
                    (GetVTrgt(s) - GetV(s)) -
                    0.5 * Power2(GetA(s)) / GetJMax(s) -
                    0.5 * Power2(GetAMin(s)) * GetInvJMaxSubInvJMin(s);
  PolynomialRoots roots = SafeCalculatePolynomialRoots(coefficients);

  max = GetA(s);
  if (roots[0] > GetA(s)) {
    max = roots[0];
  } else {
    if (roots[1] > GetA(s)) {
      max = roots[0];
    }
  }

  // The peak acceleration required to achieve the desired execution time
  // without the vtrgt constraint using PosTri and a(=0)->amin->0
  double postri_peak_acceleration_without_vtrgt =
      (desired_execution_time + GetAMin(s) * GetInvJMaxSubInvJMin(s) +
       GetA(s) / GetJMax(s)) /
      GetInvJMaxSubInvJMin(s);
  max = std::min(max, postri_peak_acceleration_without_vtrgt);

  return ClipAndVerifyRange({min, max}, {GetA(s), GetAMax(s)});
}

Range GetPosTriPeakAccelerationLimitsForPosTriZeroNegTri(
    const SyncDOFsMotionState& s, const double desired_execution_time) {
  double min = 0., max = 0.;
  // The maximum peak acceleration value can be determined by
  // calculating the peak acceleration of PosTriNegTri such
  // that t=tsync and v=vtrgt constant value

  double temp_postri_peak_acceleration_max = 0.0;
  double temp_negtri_peak_acceleration_max = 0.0;
  if (desired_execution_time > 0) {
    temp_postri_peak_acceleration_max =
        0.5 * (1 / GetInvJMaxSubInvJMin(s)) *
            (desired_execution_time + GetA(s) / GetJMax(s)) +
        ((GetVTrgt(s) - GetV(s) + 0.5 * Power2(GetA(s)) / GetJMax(s)) /
         (desired_execution_time + GetA(s) / GetJMax(s)));
    temp_negtri_peak_acceleration_max =
        temp_postri_peak_acceleration_max -
        (desired_execution_time + GetA(s) / GetJMax(s)) /
            GetInvJMaxSubInvJMin(s);

    // The peak acceleration required to achieve the desired execution time
    // without the vtrgt constraint using only PosTri
    double temp_postri_peak_acceleration_without_vtrgt =
        (desired_execution_time + GetA(s) / GetJMax(s)) /
        GetInvJMaxSubInvJMin(s);
    temp_postri_peak_acceleration_max =
        std::min(temp_postri_peak_acceleration_max,
                 temp_postri_peak_acceleration_without_vtrgt);
    temp_postri_peak_acceleration_max =
        std::min(temp_postri_peak_acceleration_max, GetAMax(s));

    // The peak acceleration required to achieve the desired execution time
    // without the vtrgt constraint using a->0 and NegTri
    double temp_negtri_peak_acceleration_without_vtrgt =
        std::min(0.0, -(desired_execution_time - GetA(s) / GetJMax(s)) /
                          GetInvJMaxSubInvJMin(s));

    temp_negtri_peak_acceleration_max =
        ClipToRange(temp_negtri_peak_acceleration_max,
                    temp_negtri_peak_acceleration_without_vtrgt, 0.0);
  }

  max = temp_postri_peak_acceleration_max;

  // a->0
  double end_velocity_a2zero = GetV(s) + 0.5 * Power2(GetA(s)) / (-GetJMin(s));

  if (end_velocity_a2zero > GetVTrgt(s)) {
    min = GetA(s);
  } else {
    // Profile: a -> apeak -> 0, so that v = vtrgt
    min = GetSqrt(
        (2.0 * (GetVTrgt(s) - GetV(s)) + (Power2(GetA(s)) / GetJMax(s))) *
        (-GetJMin(s)) * GetJMax(s) / ((-GetJMin(s)) + GetJMax(s)));
  }

  return ClipAndVerifyRange({min, max}, {GetA(s), GetAMax(s)});
}

Range GetPosTrapHoldTimeLimitsForPosTrapZeroPosTrap(
    const double desired_execution_time, const SyncDOFsMotionState& s) {
  // The maximum hold time value can be determined by calculating the hold time
  // of PosTrap and a(=0)->amax->0 such that v=vtrgt
  double max =
      (GetVTrgt(s) - GetV(s) - Power2(GetAMax(s)) * GetInvJMaxSubInvJMin(s) +
       0.5 * Power2(GetA(s)) / GetJMax(s)) /
      GetAMax(s);

  // The hold time required to achieve the desired execution time without the
  // vtrgt constraint using PosTrap and a(=0)->amax->0
  double hold_time_postrap1_without_vtrgt =
      desired_execution_time - 2.0 * GetInvJMaxSubInvJMin(s) * GetAMax(s) +
      GetA(s) / GetJMax(s);

  max = std::min(max, hold_time_postrap1_without_vtrgt);

  return ClipAndVerifyRange({0., max}, {0., desired_execution_time});
}

Range GetPosTriPeakAccelerationLimitsForPosTrapZeroPosTri(
    const double desired_execution_time, const SyncDOFsMotionState& s) {
  double max = GetAMax(s);

  // end velocity after PosTrap and a(=0)->amax->0, such that t=tsync
  double end_velocity =
      GetV(s) + GetInvJMaxSubInvJMin(s) * Power2(GetAMax(s)) -
      0.5 * Power2(GetA(s)) / GetJMax(s) +
      GetAMax(s) *
          (desired_execution_time - 2.0 * GetInvJMaxSubInvJMin(s) * GetAMax(s) +
           GetA(s) / GetJMax(s));
  if (end_velocity < GetVTrgt(s)) {
    // PosTrap and a(=0)->amax->0, such that t=tsync and v=vtrgt

    std::array<double, 3> coefficients;
    coefficients[2] = 0.5 * GetInvJMaxSubInvJMin(s);
    coefficients[1] = -GetAMax(s) * GetInvJMaxSubInvJMin(s);
    coefficients[0] = GetAMax(s) * (desired_execution_time -
                                    0.5 * GetInvJMaxSubInvJMin(s) * GetAMax(s) +
                                    GetA(s) / GetJMax(s)) -
                      (GetVTrgt(s) - GetV(s)) -
                      0.5 * Power2(GetA(s)) / GetJMax(s);

    PolynomialRoots roots = SafeCalculatePolynomialRoots(coefficients);
    max = roots[1];

    // The peak acceleration required to achieve the desired execution time
    // without the vtrgt constraint using a->amax->0 and PosTri
    double postri_peak_acceleration_without_vtrgt =
        (desired_execution_time - GetAMax(s) * GetInvJMaxSubInvJMin(s) +
         GetA(s) / GetJMax(s)) /
        GetInvJMaxSubInvJMin(s);

    max = std::min(max, postri_peak_acceleration_without_vtrgt);
    max = std::min(max, GetAMax(s));
  }

  return ClipAndVerifyRange({0., max}, {0.0, GetAMax(s)});
}

Range GetPosTriPeakAccelerationLimitsForPosTriZeroPosTrap(
    const double desired_execution_time, const SyncDOFsMotionState& s) {
  double max;
  // PosTriZeroPosTrap has two limitation possibilities for apeak1Max, but
  // does not have an limitations for apeak1Min. Hence, we have to calculate the
  // two possible limitations for apeak1Max (apeak1Max_1 and apeak1Max_2), and
  // finally we can determine the real value for apeak1Max.
  double temp1_postri_peak_acceleration_max = GetA(s);
  double temp2_postri_peak_acceleration_max = GetA(s);

  // First, we check, if the velocity difference value for the case
  // that a -> +amax -> 0 -> +amax -> 0 is greater than the desired velocity
  // difference value, vtrgt - GetV(s). If this is the case, we have to
  // calculate an alternative value for apeak1Max (apeak1Max_1) by using the
  // profile a -> +apeak -> 0 -> +amax -> 0.

  // velocity after a->amax->0->amax->0
  double end_velocity = GetV(s) + GetInvJMaxSubInvJMin(s) * Power2(GetAMax(s)) -
                        0.5 * Power2(GetA(s)) / GetJMax(s);

  // Second, we have to distinguish between two cases:
  // 1. apeak1Max = amax can be used to reach vtrgt
  // 2. apeak1Max = amax cannot be used to reach vtrgt
  double delta_time = 0.0;
  if (end_velocity <= GetVTrgt(s)) {
    temp1_postri_peak_acceleration_max = GetAMax(s);

    // We calculate the time, which we would need for the profile
    // a -> +amax -> 0 -> +amax -> hold -> 0, such that v = vtrgt
    // we would bring a to zero, then a to amax, hold, and finally back to zero,
    // such that the time for this change in velocity equals the desired time
    // difference to tsync. If this difference is lesser than the desired
    // velocity difference,

    // time for hold phase of PosTrap such that v=vtrgt after a->amax->0 and
    // PosTrap
    double hold_time_postrap =
        (GetVTrgt(s) - GetV(s) - GetInvJMaxSubInvJMin(s) * Power2(GetAMax(s)) +
         0.5 * Power2(GetA(s)) / GetJMax(s)) /
        GetAMax(s);

    // time for a->amax->0 and PosTrap such that v=vtrgt
    delta_time = 2.0 * GetInvJMaxSubInvJMin(s) * GetAMax(s) -
                 GetA(s) / GetJMax(s) + hold_time_postrap;

  } else {
    // peak acceleration in PosTri and a(=0)->amax->0 to achieve v=vtrgt
    temp1_postri_peak_acceleration_max = GetSqrt(
        (2.0 / GetInvJMaxSubInvJMin(s)) *
            (GetVTrgt(s) - GetV(s) + 0.5 * Power2(GetA(s)) / GetJMax(s)) -
        Power2(GetAMax(s)));

    // PosTri and a(=0)->amax->0 such that v=vtrgt
    // peak acceleration of PosTri needed
    double peak_acceleration = GetSqrt(
        (2.0 / GetInvJMaxSubInvJMin(s)) *
            (GetVTrgt(s) - GetV(s) + 0.5 * Power2(GetA(s)) / GetJMax(s)) -
        Power2(GetAMax(s)));

    // time for PosTri and a(=0)->amax->0 such that v=vtrgt
    delta_time = GetInvJMaxSubInvJMin(s) * (GetAMax(s) + peak_acceleration) -
                 GetA(s) / GetJMax(s);
  }

  if (delta_time > desired_execution_time) {
    // PosTri and a(=0)->amax->0, such that t=tsync and v=vtrgt
    std::array<double, 3> coefficients;
    coefficients[2] = 0.5 * GetInvJMaxSubInvJMin(s);
    coefficients[1] = -GetAMax(s) * GetInvJMaxSubInvJMin(s);
    coefficients[0] = GetAMax(s) * (desired_execution_time -
                                    0.5 * GetInvJMaxSubInvJMin(s) * GetAMax(s) +
                                    GetA(s) / GetJMax(s)) -
                      (GetVTrgt(s) - GetV(s)) -
                      0.5 * Power2(GetA(s)) / GetJMax(s);
    PolynomialRoots roots = SafeCalculatePolynomialRoots(coefficients);
    temp2_postri_peak_acceleration_max = roots[1];

    temp2_postri_peak_acceleration_max =
        ClipToRange(temp2_postri_peak_acceleration_max, {GetA(s), GetAMax(s)});
  } else {
    temp2_postri_peak_acceleration_max = GetAMax(s);
  }

  // The peak acceleration required to achieve the desired execution time
  // without the vtrgt constraint using PosTri and a->amax->0
  double postri_peak_acceleration_without_vtrgt =
      (desired_execution_time - GetAMax(s) * GetInvJMaxSubInvJMin(s) +
       GetA(s) / GetJMax(s)) /
      GetInvJMaxSubInvJMin(s);

  max = std::min(temp1_postri_peak_acceleration_max, GetAMax(s));
  max = std::min(max, temp2_postri_peak_acceleration_max);
  max = std::min(max, postri_peak_acceleration_without_vtrgt);

  return ClipAndVerifyRange({GetA(s), max}, {GetA(s), GetAMax(s)});
}

std::pair<double, bool> GetPositionErrorForPosTriZeroPosTri(
    const double desired_execution_time, const SyncDOFsMotionState& s,
    const double peak_acceleration_postri1) {
  // time durations of acceleration and deceleration phases of PosTri1
  double acceleration_time_postri1 =
      (peak_acceleration_postri1 - GetA(s)) * GetInvJMax(s);
  double deceleration_time_postri1 =
      peak_acceleration_postri1 * GetNegInvJMin(s);

  // time duration of PosTri1
  double delta_time_postri1 =
      acceleration_time_postri1 + deceleration_time_postri1;

  // State after the first postri.
  auto s1 = AUpToAPrime(peak_acceleration_postri1, s) | ADownToZero;

  // peak_acceleration_postri2 needed for PosTri2 such that v=vtrgt after
  // PosTri2
  double peak_acceleration_postri2 = 0.0;
  if (GetV(s1) < GetVTrgt(s)) {
    peak_acceleration_postri2 =
        GetSqrt(2.0 * (GetVTrgt(s) - GetV(s1)) * GetJerkConstant(s));
  }

  // time durations of acceleration and deceleration phases of PosTri2
  double acceleration_time_postri2 = peak_acceleration_postri2 * GetInvJMax(s);
  double deceleration_time_postri2 =
      peak_acceleration_postri2 * GetNegInvJMin(s);

  // time duration of PosTri
  double delta_time_postri2 =
      deceleration_time_postri2 + acceleration_time_postri2;

  // Zero
  double delta_time_zero =
      desired_execution_time - delta_time_postri1 - delta_time_postri2;

  auto s2 = s1 | HoldForDeltaT(delta_time_zero) |
            AUpToAPrime(peak_acceleration_postri2) | ADownToZero;

  return std::make_pair(
      GetPError(s2), ((GetV(s1) - GetVTrgt(s)) > kAbsVelocityErrorTolerance) ||
                         (delta_time_zero < -kAbsTimeErrorTolerance));
}

Range GetPosTriPeakAccelerationLimitsForPosTriZeroPosTri(
    const double desired_execution_time, const SyncDOFsMotionState& s) {
  double min = 0., max = 0.;
  // peak acceleration of PosTri needed to achieve v=vtrgt with a single PosTri
  double peak_acceleration_single_postri =
      GetSqrt((2.0 / GetInvJMaxSubInvJMin(s)) *
              (GetVTrgt(s) - GetV(s) + 0.5 * Power2(GetA(s)) / GetJMax(s)));

  // peak acceleration of PosTri without vtrgt constraint for a single PosTri
  double peak_acceleration_max_without_vtrgt =
      (desired_execution_time + GetA(s) / GetJMax(s)) / GetInvJMaxSubInvJMin(s);

  peak_acceleration_single_postri = std::min(
      peak_acceleration_single_postri, peak_acceleration_max_without_vtrgt);

  // If this value is greater than amax, we set the maximum allowed value for
  // apeak1 to amax.
  max = std::min(peak_acceleration_single_postri, GetAMax(s));

  // In the next step we have to consider the second PosTri face. If we would
  // bring a to zero and afterwards to amax and back to zero, we have to check,
  // whether the reached velocity is greater than our desired target velocity.
  // If this is the case, we can set the minimum allowed value for apeak1 to
  // GetA(s). Otherwise we calculate the correct value considering apeak2 =
  // amax for the second positive face.

  // velocity after a->0 and a(=0)->amax->0
  double end_velocity = GetV(s) + 0.5 * Power2(GetA(s)) / GetJMax(s) +
                        0.5 * GetInvJMaxSubInvJMin(s) * Power2(GetAMax(s));
  if (end_velocity > GetVTrgt(s)) {
    min = GetA(s);
  } else {
    // peak acceleration of PosTri needed for the profile
    // PosTri and a(=0)->amax->0 such that v=vtrgt
    min = GetSqrt(
        (2.0 / GetInvJMaxSubInvJMin(s)) *
            (GetVTrgt(s) - GetV(s) + 0.5 * Power2(GetA(s)) / GetJMax(s)) -
        Power2(GetAMax(s)));
  }

  min = std::max(min, GetA(s));
  if (min < GetA(s)) {
    min = GetA(s);
  }

  // For apeak1=apeak2
  // peak acceleration needed such that PosTriPosTri achieves v=vtrgt
  double peak_acceleration_equal_postripostri =
      GetSqrt((GetVTrgt(s) - GetV(s) + 0.5 * Power2(GetA(s)) / GetJMax(s)) /
              GetInvJMaxSubInvJMin(s));
  // time duration of the above profile
  double delta_time =
      2.0 * GetInvJMaxSubInvJMin(s) * peak_acceleration_equal_postripostri -
      GetA(s) / GetJMax(s);

  if (delta_time > desired_execution_time) {
    // Solve for PosTriPosTri such that t=tsync and v=vtrgt

    std::array<double, 3> coefficients;
    coefficients[2] = GetInvJMaxSubInvJMin(s);
    coefficients[1] = -(desired_execution_time + GetA(s) / GetJMax(s));
    coefficients[0] =
        0.5 * (1 / GetInvJMaxSubInvJMin(s)) *
            Power2(desired_execution_time + GetA(s) / GetJMax(s)) -
        (GetVTrgt(s) - GetV(s)) - 0.5 * Power2(GetA(s)) / GetJMax(s);

    PolynomialRoots roots = SafeCalculatePolynomialRoots(coefficients);

    double peak_acceleration_high = std::max(roots[0], GetA(s));
    double peak_acceleration_low = std::max(roots[1], GetA(s));

    // get position errors at these different bounds
    double abs_position_error_at_peak_acceleration_min =
        fabs(GetPositionErrorForPosTriZeroPosTri(desired_execution_time, s, min)
                 .first);
    double abs_position_error_at_peak_acceleration_low =
        fabs(GetPositionErrorForPosTriZeroPosTri(desired_execution_time, s,
                                                 peak_acceleration_low)
                 .first);
    double abs_position_error_at_peak_acceleration_high =
        fabs(GetPositionErrorForPosTriZeroPosTri(desired_execution_time, s,
                                                 peak_acceleration_high)
                 .first);
    double abs_position_error_at_peak_acceleration_max =
        fabs(GetPositionErrorForPosTriZeroPosTri(desired_execution_time, s, max)
                 .first);

    if ((abs_position_error_at_peak_acceleration_min <
             abs_position_error_at_peak_acceleration_low &&
         abs_position_error_at_peak_acceleration_min <
             abs_position_error_at_peak_acceleration_high &&
         abs_position_error_at_peak_acceleration_min <
             abs_position_error_at_peak_acceleration_max) ||
        (abs_position_error_at_peak_acceleration_low <
             abs_position_error_at_peak_acceleration_min &&
         abs_position_error_at_peak_acceleration_low <
             abs_position_error_at_peak_acceleration_high &&
         abs_position_error_at_peak_acceleration_low <
             abs_position_error_at_peak_acceleration_max)) {
      max = peak_acceleration_low;
    } else {
      min = peak_acceleration_high;
    }
  }

  return ClipAndVerifyRange({min, max}, {GetA(s), GetAMax(s)});
}

Range GetHoldAccelerationLimitsForPosLinHldPosTrap(
    const double desired_execution_time, const SyncDOFsMotionState& s) {
  double max = GetAMax(s);

  // time for hold phase of a->ahld->amax->0 such that t=tsync and v=vtrgt
  double hold_time = desired_execution_time + GetA(s) / GetJMax(s) -
                     GetInvJMaxSubInvJMin(s) * GetAMax(s);
  // hold acceleration
  if (hold_time != 0.0) {
    max = (GetVTrgt(s) - GetV(s) + 0.5 * Power2(GetA(s)) / GetJMax(s) -
           0.5 * Power2(GetAMax(s)) * GetInvJMaxSubInvJMin(s)) /
          hold_time;
    if (max > GetAMax(s)) {
      max = GetAMax(s);
    }
  }

  return ClipAndVerifyRange({GetA(s), max}, {GetA(s), GetAMax(s)});
}

Range GetHoldAccelerationLimitsForNegLinHldPosTrap(
    const double desired_execution_time, const SyncDOFsMotionState& s) {
  double min = 0., max = GetA(s);

  // velocity after a->0 and a(=0)->amax->0
  double end_velocity_a2zero_a2amax2zero =
      GetV(s) + 0.5 * Power2(GetA(s)) / GetJMax(s) +
      0.5 * GetInvJMaxSubInvJMin(s) * Power2(GetAMax(s));
  if (end_velocity_a2zero_a2amax2zero > GetVTrgt(s)) {
    // profile a->ahld->amax->hold->0, such that v = vtrgt and t = tsync

    std::array<double, 3> coefficients;
    coefficients[2] = 0.5 * GetInvJMaxSubInvJMin(s);
    coefficients[1] = -GetAMax(s) * GetInvJMaxSubInvJMin(s);
    coefficients[0] =
        -GetAMax(s) * (desired_execution_time - GetA(s) / (-GetJMin(s))) +
        (GetVTrgt(s) - GetV(s)) +
        0.5 * Power2(GetAMax(s)) * GetInvJMaxSubInvJMin(s) -
        0.5 * Power2(GetA(s)) / (-GetJMin(s));

    PolynomialRoots roots = SafeCalculatePolynomialRoots(coefficients);
    min = roots[1];

  } else {
    // profile a->0->amax->hold->0, such that v = vtrgt
    double delta_time =
        GetA(s) / (-GetJMin(s)) + GetAMax(s) * GetInvJMaxSubInvJMin(s) +
        (GetVTrgt(s) - GetV(s) - 0.5 * Power2(GetA(s)) / (-GetJMin(s)) -
         0.5 * Power2(GetAMax(s)) * GetInvJMaxSubInvJMin(s)) /
            GetAMax(s);

    if (delta_time > desired_execution_time) {
      // profile a->ahld->amax->hold->0, such that v = vtrgt and t = tsync

      std::array<double, 3> coefficients;
      coefficients[2] = 0.5 * GetInvJMaxSubInvJMin(s);
      coefficients[1] = -GetAMax(s) * GetInvJMaxSubInvJMin(s);
      coefficients[0] =
          -GetAMax(s) * (desired_execution_time - GetA(s) / (-GetJMin(s))) +
          (GetVTrgt(s) - GetV(s)) +
          0.5 * Power2(GetAMax(s)) * GetInvJMaxSubInvJMin(s) -
          0.5 * Power2(GetA(s)) / (-GetJMin(s));
      PolynomialRoots roots = SafeCalculatePolynomialRoots(coefficients);
      min = roots[1];
    }
  }

  // Profile: a->ahld->amax->0, such that t = tsync and v = vtrgt
  double hold_acceleration_without_vtrgt =
      GetAMax(s) - (1 / GetInvJMaxSubInvJMin(s)) *
                       (desired_execution_time - GetA(s) / (-GetJMin(s)));

  min = std::max(min, hold_acceleration_without_vtrgt);

  // Profile: a->hold->amax->0, such that t = tsync
  // end velocity after a->hold->amax->0, such that t = tsync
  double end_velocity_ahold2amax2zero =
      GetV(s) +
      GetA(s) * (desired_execution_time - GetAMax(s) * GetInvJMaxSubInvJMin(s) +
                 GetA(s) / GetJMax(s)) +
      0.5 * Power2(GetAMax(s)) * GetInvJMaxSubInvJMin(s) -
      0.5 * Power2(GetA(s)) / GetJMax(s);
  if (end_velocity_ahold2amax2zero <= GetVTrgt(s)) {
    // Profile: a->ahld->hold->amax->0, such that t = tsync and v = vtrgt

    std::array<double, 3> coefficients;
    coefficients[2] = 0.5 * GetInvJMaxSubInvJMin(s);
    coefficients[1] = desired_execution_time -
                      GetAMax(s) * GetInvJMaxSubInvJMin(s) -
                      GetA(s) / (-GetJMin(s));
    coefficients[0] = 0.5 * GetInvJMaxSubInvJMin(s) * Power2(GetAMax(s)) +
                      0.5 * Power2(GetA(s)) / (-GetJMin(s)) -
                      (GetVTrgt(s) - GetV(s));

    PolynomialRoots roots = SafeCalculatePolynomialRoots(coefficients);
    if (roots[0] > 0.0) {
      max = roots[0];
    } else {
      if (roots[1] > 0.0) {
        max = roots[1];
      }
    }
  }

  return ClipAndVerifyRange({min, max}, {0., GetA(s)});
}

Range GetPosTrapHoldTimeLimitsForPosTrapHldNegLin(
    const double desired_execution_time, const SyncDOFsMotionState& s) {
  // Profile: a->amax->hold->0 such that v=vtrgt
  double max = (GetVTrgt(s) - GetV(s) -
                0.5 * Power2(GetAMax(s)) * GetInvJMaxSubInvJMin(s) +
                0.5 * Power2(GetA(s)) / GetJMax(s)) /
               GetAMax(s);

  // Profile: a->amax->hold->0 such that t=tsync
  double hold_time_postrap_without_vtrgt =
      desired_execution_time - GetInvJMaxSubInvJMin(s) * GetAMax(s) +
      GetA(s) / GetJMax(s);

  max = std::min(max, hold_time_postrap_without_vtrgt);

  return ClipAndVerifyRange({0., max}, {0., desired_execution_time});
}

Range GetPosTriPeakAccelerationLimitsForPosTriHldNegLin(
    const double desired_execution_time, const SyncDOFsMotionState& s) {
  // Profile: a->apeak->0 such that v = vtrgt
  double max =
      GetSqrt((2.0 / GetInvJMaxSubInvJMin(s)) *
              (GetVTrgt(s) - GetV(s) + 0.5 * Power2(GetA(s)) / GetJMax(s)));

  // Profile: a->apeak->0 such that t = tsync
  double peak_acceleration_max_without_vtrgt =
      (desired_execution_time + GetA(s) / GetJMax(s)) / GetInvJMaxSubInvJMin(s);

  max = std::min(max, peak_acceleration_max_without_vtrgt);

  // end velocity after a->hold->0 such that t=tsync
  double end_velocity = GetV(s) + GetA(s) * desired_execution_time -
                        0.5 * Power2(GetA(s)) / (-GetJMin(s));

  double min = GetA(s);
  if (end_velocity <= GetVTrgt(s)) {
    // Profile a->ahold->hold->0 such that t=tsync and v=vtrgt

    std::array<double, 3> coefficients;
    coefficients[2] = -0.5 * GetInvJMaxSubInvJMin(s);
    coefficients[1] = desired_execution_time + GetA(s) / GetJMax(s);
    coefficients[0] =
        -0.5 * Power2(GetA(s)) / GetJMax(s) - (GetVTrgt(s) - GetV(s));

    PolynomialRoots roots = SafeCalculatePolynomialRoots(coefficients);
    if ((roots[1] > GetA(s)) && (roots[1] < GetAMax(s))) {
      min = roots[1];
    } else {
      if ((roots[0] > GetA(s)) && (roots[0] < GetAMax(s))) {
        min = roots[0];
      }
    }
  }

  return ClipAndVerifyRange({min, max}, {GetA(s), GetAMax(s)});
}

Range GetHoldAccelerationLimitsForNegLinHldPosTri(
    const double desired_execution_time, const SyncDOFsMotionState& s) {
  double min = 0.;
  // Time for a->0->apeak->0 such that v = vtrgt
  double delta_time =
      GetA(s) / (-GetJMin(s)) +
      GetInvJMaxSubInvJMin(s) * GetSqrt((2.0 * (GetVTrgt(s) - GetV(s)) -
                                         Power2(GetA(s)) / (-GetJMin(s))) /
                                        GetInvJMaxSubInvJMin(s));
  if (delta_time < desired_execution_time) {
    // Profile a -> apeak -> 0 such that v = vtrgt
    double peak_acceleration =
        GetSqrt((2.0 * (GetVTrgt(s) - GetV(s)) + Power2(GetA(s)) / GetJMax(s)) /
                GetInvJMaxSubInvJMin(s));

    if (peak_acceleration > GetAMax(s)) {
      // Profile: a->ahld->hold->amax->0, such that t = tsync and v = vtrgt
      std::array<double, 3> coefficients;
      coefficients[2] = 0.5 * GetInvJMaxSubInvJMin(s);
      coefficients[1] = desired_execution_time -
                        GetAMax(s) * GetInvJMaxSubInvJMin(s) -
                        GetA(s) / (-GetJMin(s));
      coefficients[0] = 0.5 * GetInvJMaxSubInvJMin(s) * Power2(GetAMax(s)) +
                        0.5 * Power2(GetA(s)) / (-GetJMin(s)) -
                        (GetVTrgt(s) - GetV(s));
      PolynomialRoots roots = CalculatePolynomialRoots(coefficients);
      if (!roots.empty()) {
        if (roots[0] > 0.0) {
          min = roots[0];
        } else {
          if (roots[1] > 0.0) {
            min = roots[1];
          }
        }
      }
    }
  } else {
    // The following minimum ahold value is calculated by using the profile
    // a->ahld->apeak->0 such that t=tsync and v=vtrgt
    double diff_time = desired_execution_time - GetA(s) / (-GetJMin(s));
    if (diff_time != 0) {
      min =
          0.5 *
          (2.0 *
               (GetVTrgt(s) - GetV(s) - 0.5 * Power2(GetA(s)) / (-GetJMin(s))) /
               diff_time -
           diff_time / GetInvJMaxSubInvJMin(s));
    }
  }

  // end velocity after a->hold->0 such that t=tsync
  double end_velocity = GetV(s) + GetA(s) * desired_execution_time -
                        0.5 * Power2(GetA(s)) / (-GetJMin(s));

  double max = GetA(s);
  if (end_velocity > GetVTrgt(s)) {
    // Profile a->ahold->hold->0 such that t=tsync and v=vtrgt
    double hold_time = desired_execution_time - GetA(s) / (-GetJMin(s));

    if (hold_time != 0) {
      max = (GetVTrgt(s) - GetV(s) - 0.5 * Power2(GetA(s)) / (-GetJMin(s))) /
            hold_time;
    }
  }

  return ClipAndVerifyRange({min, max}, {0., GetA(s)});
}

Range GetPosTriPeakAccelerationLimitsForPosLinHldPosTri(
    const double desired_execution_time, const SyncDOFsMotionState& s) {
  // calculate apeak1 for the case a->hold->apeak1->0 such that v = vtrgt and
  // t = tsync;
  std::array<double, 3> coefficients;
  coefficients[2] = 0.5 * GetInvJMaxSubInvJMin(s);
  coefficients[1] = -GetA(s) * GetInvJMaxSubInvJMin(s);
  coefficients[0] = GetA(s) * desired_execution_time +
                    0.5 * Power2(GetA(s)) / GetJMax(s) -
                    (GetVTrgt(s) - GetV(s));
  PolynomialRoots roots = SafeCalculatePolynomialRoots(coefficients);
  double max = roots[0];

  // calculate apeak1 for only PosTri without vtrgt constraint such that t=tsync
  double peak_acceleration_max_without_vtrgt =
      (desired_execution_time + GetA(s) / GetJMax(s)) / GetInvJMaxSubInvJMin(s);

  max = std::min(max, peak_acceleration_max_without_vtrgt);

  // profile: a->ahld->hold->0 such that v=vtrgt and t=tsync

  coefficients[2] = 0.5 * GetInvJMaxSubInvJMin(s);
  coefficients[1] = -(desired_execution_time + GetA(s) / GetJMax(s));
  coefficients[0] =
      0.5 * Power2(GetA(s)) / GetJMax(s) + (GetVTrgt(s) - GetV(s));

  roots = SafeCalculatePolynomialRoots(coefficients);
  double min = GetA(s);
  if (roots[1] < GetAMax(s)) {
    min = roots[1];
  } else {
    if (roots[0] < GetAMax(s)) {
      min = roots[0];
    }
  }

  return ClipAndVerifyRange({min, max}, {GetA(s), GetAMax(s)});
}

std::pair<double, bool> GetPositionErrorForPosTrapZeroNegTrap(
    const std::array<double, 3>& position_error_poly_coeff,
    const std::array<double, 2>& delta_velocity_hold_phase_negtrap_poly_coeff,
    const std::array<double, 2>& zero_time_poly_coeff,
    const SyncDOFsMotionState& s, const double hold_time_postrap) {
  // check if solution is invalid
  bool invalid_solution =
      ((EvalPolynomial(delta_velocity_hold_phase_negtrap_poly_coeff,
                       hold_time_postrap) > kAbsVelocityErrorTolerance) ||
       (EvalPolynomial(zero_time_poly_coeff, hold_time_postrap) <
        -kAbsTimeErrorTolerance));

  // return position error
  return std::make_pair(
      EvalPolynomial(position_error_poly_coeff, hold_time_postrap),
      invalid_solution);
}

std::pair<double, bool> GetPositionErrorForPosTrapZeroNegTri(
    const std::array<double, 5>& position_error_poly_coeff,
    const std::array<double, 3>& delta_velocity_hold_phase_postrap_poly_coeff,
    const std::array<double, 3>& zero_time_poly_coeff,
    const SyncDOFsMotionState& s, const double peak_acceleration_negtri) {
  // check if solution is invalid
  bool invalid_solution =
      ((EvalPolynomial(delta_velocity_hold_phase_postrap_poly_coeff,
                       peak_acceleration_negtri) <
        -kAbsVelocityErrorTolerance) ||
       (EvalPolynomial(zero_time_poly_coeff, peak_acceleration_negtri) <
        -kAbsTimeErrorTolerance));

  // return position error
  return std::make_pair(
      EvalPolynomial(position_error_poly_coeff, peak_acceleration_negtri),
      invalid_solution);
}

std::pair<double, bool> GetPositionErrorForPosTriZeroNegTrap(
    const std::array<double, 5>& position_error_poly_coeff,
    const std::array<double, 3>& delta_velocity_hold_phase_negtrap_poly_coeff,
    const std::array<double, 3>& zero_time_poly_coeff,
    const SyncDOFsMotionState& s, const double peak_acceleration_postri) {
  // check if solution is invalid
  bool invalid_solution =
      ((EvalPolynomial(delta_velocity_hold_phase_negtrap_poly_coeff,
                       peak_acceleration_postri) >
        kAbsVelocityErrorTolerance) ||
       (EvalPolynomial(zero_time_poly_coeff, peak_acceleration_postri) <
        -kAbsTimeErrorTolerance));

  // return position error
  return std::make_pair(
      EvalPolynomial(position_error_poly_coeff, peak_acceleration_postri),
      invalid_solution);
}

std::pair<double, bool> GetPositionErrorForPosTriZeroNegTri(
    const SyncDOFsMotionState& s, const double peak_acceleration_postri) {
  auto s1 = AUpToAPrimePosTri(peak_acceleration_postri, s);

  double delta_time_zero = GetTSync(s) - GetT(VDownToVTrgtNegTri(s1));

  return std::make_pair(
      HoldForDeltaT(delta_time_zero, s1) | VDownToVTrgtNegTri | GetPError,
      ((GetV(s1) - GetVTrgt(s)) < -kAbsVelocityErrorTolerance) ||
          (delta_time_zero < -kAbsTimeErrorTolerance));
}

std::pair<double, bool> GetPositionErrorForPosTrapZeroPosTrap(
    const std::array<double, 2>& position_error_poly_coeff,
    const std::array<double, 2>& delta_velocity_hold_phase_postrap2_poly_coeff,
    const SyncDOFsMotionState& s, const double hold_time_postrap1) {
  // check if solution is invalid
  bool invalid_solution =
      (EvalPolynomial(delta_velocity_hold_phase_postrap2_poly_coeff,
                      hold_time_postrap1) < -kAbsVelocityErrorTolerance);

  // return position error
  return std::make_pair(
      EvalPolynomial(position_error_poly_coeff, hold_time_postrap1),
      invalid_solution);
}

std::pair<double, bool> GetPositionErrorForPosTrapZeroPosTri(
    const std::array<double, 5>& position_error_poly_coeff,
    const std::array<double, 3>& delta_velocity_hold_phase_postrap_poly_coeff,
    const std::array<double, 3>& zero_time_poly_coeff,
    const SyncDOFsMotionState& s, const double peak_acceleration_postri) {
  // check if solution is invalid
  bool invalid_solution =
      ((EvalPolynomial(delta_velocity_hold_phase_postrap_poly_coeff,
                       peak_acceleration_postri) <
        -kAbsVelocityErrorTolerance) ||
       (EvalPolynomial(zero_time_poly_coeff, peak_acceleration_postri) <
        -kAbsTimeErrorTolerance));

  // return position error
  return std::make_pair(
      EvalPolynomial(position_error_poly_coeff, peak_acceleration_postri),
      invalid_solution);
}

std::pair<double, bool> GetPositionErrorForPosTriZeroPosTrap(
    const std::array<double, 5>& position_error_poly_coeff,
    const std::array<double, 3>& delta_velocity_hold_phase_postrap_poly_coeff,
    const std::array<double, 3>& zero_time_poly_coeff,
    const SyncDOFsMotionState& s, const double peak_acceleration_postri) {
  // check if solution is invalid
  bool invalid_solution =
      (EvalPolynomial(delta_velocity_hold_phase_postrap_poly_coeff,
                      peak_acceleration_postri) < -kAbsVelocityErrorTolerance ||
       EvalPolynomial(zero_time_poly_coeff, peak_acceleration_postri) <
           -kAbsTimeErrorTolerance);

  // return position error
  return std::make_pair(
      EvalPolynomial(position_error_poly_coeff, peak_acceleration_postri),
      invalid_solution);
}

std::pair<double, bool> GetPositionErrorForPosLinHldPosTrap(
    const double desired_execution_time, const SyncDOFsMotionState& s,
    const double hold_acceleration) {
  // time duration of PosLin
  double delta_time_poslin = (hold_acceleration - GetA(s)) * GetInvJMax(s);

  // Hld
  double hold_time_ahld = 0.0;
  if (hold_acceleration != GetAMax(s)) {
    hold_time_ahld =
        (GetAMax(s) * desired_execution_time +
         GetA(s) * GetAMax(s) * GetInvJMax(s) - (GetVTrgt(s) - GetV(s)) -
         0.5 * Power2(GetAMax(s)) * GetInvJerkConstant(s) -
         0.5 * Power2(GetA(s)) * GetInvJMax(s)) /
        (GetAMax(s) - hold_acceleration);
  }

  // PosTrap

  // acceleration time of Postrap
  double acceleration_time_postrap =
      (GetAMax(s) - hold_acceleration) * GetInvJMax(s);

  // hold time of PosTrap
  double hold_time_postrap = desired_execution_time - delta_time_poslin -
                             hold_time_ahld - acceleration_time_postrap -
                             GetAMaxToZeroDeltaT(s);

  return std::make_pair(AUpToAPrime(hold_acceleration, s) |
                            HoldForDeltaT(hold_time_ahld) | VUpToVTrgtPosTrap |
                            GetPError,
                        (hold_time_ahld < -kAbsTimeErrorTolerance) ||
                            (hold_time_postrap < -kAbsTimeErrorTolerance));
}

std::pair<double, bool> GetPositionErrorForNegLinHldPosTrap(
    const double desired_execution_time, const SyncDOFsMotionState& s,
    const double hold_acceleration) {
  // Hld
  double hold_time_ahld = 0.0;
  if (hold_acceleration != GetAMax(s)) {
    hold_time_ahld = (GetAMax(s) * (desired_execution_time -
                                    GetInvJerkConstant(s) *
                                        (GetAMax(s) - hold_acceleration) -
                                    GetA(s) / (-GetJMin(s))) -
                      (GetVTrgt(s) - GetV(s)) +
                      0.5 * GetInvJerkConstant(s) *
                          (Power2(GetAMax(s)) - Power2(hold_acceleration)) +
                      0.5 * Power2(GetA(s)) * GetNegInvJMin(s)) /
                     (GetAMax(s) - hold_acceleration);
  }

  // acceleration and deceleration times of Postrap
  double acceleration_time_postrap =
      (GetAMax(s) - hold_acceleration) * GetInvJMax(s);

  // time duration of NegLin
  double delta_time_neglin = (GetA(s) - hold_acceleration) * GetNegInvJMin(s);

  // hold time of PosTrap
  double hold_time_postrap = desired_execution_time - delta_time_neglin -
                             hold_time_ahld - acceleration_time_postrap -
                             GetAMaxToZeroDeltaT(s);

  return std::make_pair(ADownToAPrime(hold_acceleration, s) |
                            HoldForDeltaT(hold_time_ahld) | AUpToAMax |
                            HoldToTSyncAndThen(AMaxDownToZero) | GetPError,
                        (hold_time_ahld < -kAbsTimeErrorTolerance) ||
                            (hold_time_postrap < -kAbsTimeErrorTolerance));
}

std::pair<double, bool> GetPositionErrorForPosTrapHldNegLin(
    const double desired_execution_time, const SyncDOFsMotionState& s,
    const double hold_time_postrap) {
  // Hld
  double hold_time_ahld = desired_execution_time - hold_time_postrap -
                          GetInvJerkConstant(s) * GetAMax(s) +
                          GetA(s) * GetInvJMax(s);

  return std::make_pair(
      AUpToAMax(s) | HoldForDeltaT(hold_time_postrap) |
          ADownToAHldThenHoldToVTrgtTsyncAndThen(ADownToZero) | GetPError,
      (hold_time_ahld < -kAbsTimeErrorTolerance));
}

std::pair<double, bool> GetPositionErrorForPosTriHldNegLin(
    const double desired_execution_time, const SyncDOFsMotionState& s,
    const double peak_acceleration_postri) {
  // Hld
  double hold_time_ahld = desired_execution_time -
                          GetInvJerkConstant(s) * peak_acceleration_postri +
                          GetA(s) * GetInvJMax(s);

  return std::make_pair(
      AUpToAPrime(peak_acceleration_postri, s) |
          ADownToAHldThenHoldToVTrgtTsyncAndThen(ADownToZero) | GetPError,
      (hold_time_ahld < -kAbsTimeErrorTolerance));
}

std::pair<double, bool> GetPositionErrorForPosLinHldPosTri(
    const double desired_execution_time, const SyncDOFsMotionState& s,
    const double peak_acceleration_postri) {
  // hold time
  double hold_time_ahld = desired_execution_time -
                          peak_acceleration_postri * GetInvJerkConstant(s) +
                          GetA(s) * GetInvJMax(s);

  if (hold_time_ahld < kAbsZeroTimeEpsilon) {
    hold_time_ahld = 0.0;
  }

  // hold acceleration
  double hold_acceleration = GetA(s);
  if (hold_time_ahld > 0) {
    hold_acceleration =
        (GetVTrgt(s) - GetV(s) -
         0.5 * GetInvJerkConstant(s) * Power2(peak_acceleration_postri) +
         0.5 * Power2(GetA(s)) * GetInvJMax(s)) /
        hold_time_ahld;
  }

  return std::make_pair(
      AUpToAPrime(hold_acceleration, s) | HoldForDeltaT(hold_time_ahld) |
          AUpToAPrimePosTri(peak_acceleration_postri) | GetPError,
      (hold_time_ahld < -kAbsTimeErrorTolerance) ||
          ((hold_acceleration - peak_acceleration_postri) >
           kAbsAccelerationErrorTolerance));
}

std::pair<double, bool> GetPositionErrorForNegLinHldPosTri(
    const double desired_execution_time, const SyncDOFsMotionState& s,
    const double hold_acceleration) {
  // peak acceleration

  double coeff_a = 0.5 * GetInvJerkConstant(s);
  double coeff_b = -hold_acceleration * GetInvJerkConstant(s);
  double coeff_c =
      hold_acceleration * (desired_execution_time +
                           0.5 * GetInvJerkConstant(s) * hold_acceleration -
                           GetA(s) * GetNegInvJMin(s)) +
      0.5 * Power2(GetA(s)) * GetNegInvJMin(s) - (GetVTrgt(s) - GetV(s));
  double inv_twice_coeff_a = GetJerkConstant(s);

  double peak_acceleration_postri =
      (-coeff_b + GetSqrt(Power2(coeff_b) - 4 * coeff_a * coeff_c)) *
      inv_twice_coeff_a;

  // Hld
  double hold_time_ahld =
      desired_execution_time -
      GetInvJerkConstant(s) * (peak_acceleration_postri - hold_acceleration) -
      GetA(s) * GetNegInvJMin(s);

  return std::make_pair(
      ADownToAPrime(hold_acceleration, s) | HoldForDeltaT(hold_time_ahld) |
          AUpToAPrimePosTri(peak_acceleration_postri) | GetPError,
      (hold_time_ahld < -kAbsTimeErrorTolerance));
}

// Template Anderson-Bjorck-King (ABK) Method to Find Roots over F within the
// given limits.
template <typename F>
double GetRootABKMethod(const F& error_function, const Range& limits) {
  double lower_limit = limits.min;
  double upper_limit = limits.max;

  auto [position_error_at_lower_limit, invalid_solution_at_lower_limit] =
      error_function(lower_limit);

  if ((fabs(position_error_at_lower_limit) < kABKEpsilon) &&
      !invalid_solution_at_lower_limit) {
    return lower_limit;
  }

  auto [position_error_at_upper_limit, invalid_solution_at_upper_limit] =
      error_function(upper_limit);

  if ((fabs(position_error_at_upper_limit) < kABKEpsilon) &&
      !invalid_solution_at_upper_limit) {
    return upper_limit;
  }

  if (GetSign(position_error_at_lower_limit) ==
      GetSign(position_error_at_upper_limit)) {
    if (fabs(position_error_at_lower_limit) <
        fabs(position_error_at_upper_limit)) {
      if (!invalid_solution_at_lower_limit) {
        return lower_limit;
      }
    } else {
      if (!invalid_solution_at_upper_limit) {
        return upper_limit;
      }
    }
  }

  double param_1 = lower_limit;
  double param_2 = upper_limit;

  double position_error_at_param_1 = position_error_at_lower_limit;
  double position_error_at_param_2 = position_error_at_upper_limit;

  for (unsigned int x = 0; x < kNumberOfLoopsABKBlue; x++) {
    if (fabs(param_2 - param_1) <= kABKEpsilon) {
      break;
    }

    // bisection step
    double param_3 = 0.5 * (param_1 + param_2);
    auto [position_error_at_param_3, invalid_solution_at_param_3] =
        error_function(param_3);

    if ((fabs(position_error_at_param_3) < kABKEpsilon) &&
        !invalid_solution_at_param_3) {
      return param_3;
    }

    // interval determination step 1
    if ((position_error_at_param_3 * position_error_at_param_2) < 0.0) {
      param_1 = param_2;
      param_2 = param_3;
      position_error_at_param_1 = position_error_at_param_2;
      position_error_at_param_2 = position_error_at_param_3;
    } else {
      param_2 = param_3;
      position_error_at_param_2 = position_error_at_param_3;
    }

    // secant step
    double secant_12 = (position_error_at_param_1 - position_error_at_param_2) /
                       (param_1 - param_2);
    // check for valid secant value
    // ignore the secant step if not valid
    // this SHOULD NOT happen since  the function is monotically increasing or
    // decreasing within the limits
    // but still checking here for completeness
    if (fabs(secant_12) <= kABKFunctionEpsilon) {
      continue;
    }
    param_3 = param_2 - position_error_at_param_2 / secant_12;
    std::tie(position_error_at_param_3, invalid_solution_at_param_3) =
        error_function(param_3);

    if ((fabs(position_error_at_param_3) < kABKEpsilon) &&
        !invalid_solution_at_param_3) {
      return param_3;
    }

    // determination of new inclusion interval
    if ((position_error_at_param_3 * position_error_at_param_2) < 0.0) {
      param_1 = param_2;
      param_2 = param_3;
      position_error_at_param_1 = position_error_at_param_2;
      position_error_at_param_2 = position_error_at_param_3;
    } else {
      double g = 1.0 - (position_error_at_param_3 / position_error_at_param_2);
      if (g <= 0.0) {
        g = 0.5;
      }
      param_2 = param_3;
      position_error_at_param_1 *= g;
      position_error_at_param_2 = position_error_at_param_3;
    }
  }

  double param_root =
      ((fabs(position_error_at_param_2) < fabs(position_error_at_param_1))
           ? (param_2)
           : (param_1));
  auto [position_error_at_param_root, invalid_solution_at_param_root] =
      error_function(param_root);

  // check if root is invalid
  if (invalid_solution_at_param_root) {
    // check if the limits are valid
    if (!invalid_solution_at_lower_limit && !invalid_solution_at_upper_limit) {
      // choose the limit with the smallest position error
      if (fabs(position_error_at_lower_limit) <
          fabs(position_error_at_upper_limit)) {
        return lower_limit;
      } else {
        return upper_limit;
      }
    } else {
      // choose the valid limit, if any
      if (!invalid_solution_at_lower_limit) {
        return lower_limit;
      } else {
        if (!invalid_solution_at_upper_limit) {
          return upper_limit;
        }
      }
    }

  } else {
    // if the root is valid, still check with the position errors of valid
    // limits
    if (!invalid_solution_at_lower_limit &&
        (fabs(position_error_at_lower_limit) <
         fabs(position_error_at_param_root))) {
      return lower_limit;
    }
    if (!invalid_solution_at_upper_limit &&
        (fabs(position_error_at_upper_limit) <
         fabs(position_error_at_param_root))) {
      return upper_limit;
    }
  }

  return param_root;
}

constexpr auto FinishBlueProfile = Curry<SyncDOFsMotionState, Profile,
                                         IsMotionState>(
    [](const SyncDOFsMotionState& s_orig, Profile profile, const auto& s) {
      // Check for failure
      if ((fabs(GetP(s) - GetPTrgt(s)) >
           (fabs(kRelPositionErrorTolerance * (GetP(s_orig) - GetPTrgt(s))) +
            kAbsPositionErrorTolerance)) ||
          (fabs(GetV(s) - GetVTrgt(s)) >
           (fabs(kRelVelocityErrorTolerance * (GetV(s_orig) - GetVTrgt(s))) +
            kAbsVelocityErrorTolerance))) {
        if (GetTSync(s) < kMaxAllowedSyncTimeSeconds) {
          return SetFailure(s_orig);
        }
      }

      // Finish off the polynomials by leaving one that goes on forever
      double flip_mult = IsFlipped(s) ? -1.0 : 1.0;
      GetPolynomials(s)->TrimSegmentCount(GetPolynomialIndex(s));
      GetPolynomials(s)->AddSegment({.position = flip_mult * GetPTrgt(s),
                                     .velocity = flip_mult * GetVTrgt(s),
                                     .acceleration = 0.,
                                     .jerk = 0.,
                                     .start_time = GetT(s),
                                     .end_time = kInfinity});
      return SetPolynomialIndex(GetPolynomials(s)->GetSegmentCount(), s) |
             SetProfile(profile);
    });

}  // namespace

SyncDOFsMotionState CalcTrajPosTrapZeroNegTrap(const SyncDOFsMotionState& s) {
  double desired_execution_time = GetTSync(s) - GetT(s);
  // compute polynomial coefficients

  double amax = GetAMax(s);
  double jmax = GetJMax(s);

  double amin = GetAMin(s);
  double jmin = GetJMin(s);

  double vtrgt = GetVTrgt(s);
  double ptrgt = GetPTrgt(s);

  double texe = desired_execution_time;

  // postrap hold time
  std::array<double, 3> position_error_poly_coeff;
  position_error_poly_coeff[2] = amax * (-0.5 + (0.5 * amax) / amin);
  position_error_poly_coeff[1] =
      (amax *
       (amax * amin * (0.5 * jmax - jmin) +
        Power2(amax) * (-0.5 * jmax + 0.5 * jmin) +
        jmin * (-0.5 * Power2(GetA(s)) + GetA(s) * amin + 0.5 * Power2(amin) +
                amin * jmax * texe + jmax * GetV(s) - jmax * vtrgt))) /
      (amin * jmax * jmin);
  position_error_poly_coeff[0] =
      (Power3(amax) * amin *
           (-kOneSixth * Power2(jmax) + 0.5 * jmax * jmin -
            kOneThird * Power2(jmin)) +
       Power4(amin) *
           (kOneTwentyFourth * Power2(jmax) - kOneTwentyFourth * Power2(jmin)) +
       Power4(amax) *
           (0.125 * Power2(jmax) - 0.25 * jmax * jmin + 0.125 * Power2(jmin)) +
       amin * Power2(jmin) *
           (-kOneSixth * Power3(GetA(s)) - 0.5 * Power2(GetA(s)) * jmax * texe +
            Power2(jmax) * (GetP(s) - ptrgt + texe * GetV(s))) +
       Power2(amin) * Power2(jmin) *
           (-0.25 * Power2(GetA(s)) + 0.5 * jmax * GetV(s) -
            0.5 * jmax * vtrgt) +
       Power2(amax) * jmin *
           (Power2(GetA(s)) * (0.25 * jmax - 0.25 * jmin) +
            Power2(amin) * (-0.25 * jmax + 0.25 * jmin) +
            GetA(s) * amin * (-0.5 * jmax + 0.5 * jmin) +
            amin * jmax * (-0.5 * jmax + 0.5 * jmin) * texe +
            jmax * (-0.5 * jmax * GetV(s) + 0.5 * jmin * GetV(s) +
                    0.5 * jmax * vtrgt - 0.5 * jmin * vtrgt)) +
       Power2(jmin) * (0.125 * Power4(GetA(s)) +
                       Power2(GetA(s)) * jmax * (-0.5 * GetV(s) + 0.5 * vtrgt) +
                       Power2(jmax) * (0.5 * Power2(GetV(s)) - GetV(s) * vtrgt +
                                       0.5 * Power2(vtrgt)))) /
      (amin * Power2(jmax) * Power2(jmin));

  // negtrap hold phase delta velocity
  // linear polynomial coefficients
  // dv6 = a*x + b
  std::array<double, 2> delta_velocity_hold_phase_negtrap_poly_coeff;
  delta_velocity_hold_phase_negtrap_poly_coeff[1] = -amax;
  delta_velocity_hold_phase_negtrap_poly_coeff[0] =
      (0.5 *
       (Power2(amax) * (jmax - jmin) + Power2(amin) * (-jmax + jmin) +
        jmin * (Power2(GetA(s)) - 2. * jmax * GetV(s) + 2. * jmax * vtrgt))) /
      (jmax * jmin);

  // zero time
  // linear polynomial coefficients
  // dt4 = a*x + b
  std::array<double, 2> zero_time_poly_coeff;
  zero_time_poly_coeff[1] = -1 + amax / amin;
  zero_time_poly_coeff[0] =
      (0.5 *
       (2. * amax * amin * (jmax - jmin) + Power2(amax) * (-jmax + jmin) +
        Power2(amin) * (-jmax + jmin) +
        2. * amin * jmin * (GetA(s) + jmax * texe) -
        jmin * (Power2(GetA(s)) - 2. * jmax * GetV(s) + 2. * jmax * vtrgt))) /
      (amin * jmax * jmin);

  // Find solution using bisection method
  double hold_time_postrap = GetRootABKMethod(
      absl::bind_front(&GetPositionErrorForPosTrapZeroNegTrap,
                       std::cref(position_error_poly_coeff),
                       std::cref(delta_velocity_hold_phase_negtrap_poly_coeff),
                       std::cref(zero_time_poly_coeff), std::cref(s)),
      GetPosTrapHoldTimeLimitsForPosTrapZeroNegTrap(desired_execution_time, s));

  auto [position_error_optimal, invalid_solution] =
      GetPositionErrorForPosTrapZeroNegTrap(
          position_error_poly_coeff,
          delta_velocity_hold_phase_negtrap_poly_coeff, zero_time_poly_coeff, s,
          hold_time_postrap);
  if (invalid_solution ||
      (fabs(position_error_optimal) >
       (fabs(kRelPositionErrorTolerance * (GetP(s) - GetPTrgt(s))) +
        kAbsPositionErrorTolerance))) {
    if ((invalid_solution) || (GetTSync(s) < kMaxAllowedSyncTimeSeconds)) {
      return SetFailure(s);
    }
  }
  hold_time_postrap = std::max(hold_time_postrap, 0.);

  return GeneratePolynomials(AUpToAMax | HoldForDeltaT(hold_time_postrap) |
                                 AMaxDownToZero |
                                 HoldToTSyncAndThen(VDownToVTrgtNegTrap),
                             s) |
         FinishBlueProfile(s, Profile::kPosTrapZeroNegTrap);
}

SyncDOFsMotionState CalcTrajPosTrapZeroNegTri(const SyncDOFsMotionState& s) {
  double desired_execution_time = GetTSync(s) - GetT(s);
  // compute polynomial coefficients

  double amax = GetAMax(s);
  double jmax = GetJMax(s);

  double jmin = GetJMin(s);

  double vtrgt = GetVTrgt(s);
  double ptrgt = GetPTrgt(s);

  double texe = desired_execution_time;

  // negtri peak acceleration
  // quartic polynomial coefficients
  // dp = a*x^4 + b*x^3 + c*x^2 + d*x + e
  std::array<double, 5> position_error_poly_coeff;
  position_error_poly_coeff[4] =
      (-0.125 * Power2(jmax) + 0.25 * jmax * jmin - 0.125 * Power2(jmin)) /
      (amax * Power2(jmax) * Power2(jmin));
  position_error_poly_coeff[3] =
      kOneThird / Power2(jmax) + kOneSixth / Power2(jmin) - 0.5 / (jmax * jmin);
  position_error_poly_coeff[2] =
      (Power2(GetA(s)) * (0.25 * jmax - 0.25 * jmin) +
       Power2(amax) * (0.25 * jmax - 0.25 * jmin) +
       GetA(s) * amax * (-0.5 * jmax + 0.5 * jmin) +
       amax * jmax * (-0.5 * jmax + 0.5 * jmin) * texe +
       jmax * (-0.5 * jmax * GetV(s) + 0.5 * jmin * GetV(s) +
               0.5 * jmax * vtrgt - 0.5 * jmin * vtrgt)) /
      (amax * Power2(jmax) * jmin);
  position_error_poly_coeff[1] = 0.0;
  position_error_poly_coeff[0] =
      (Power4(amax) * (-kOneTwentyFourth * Power2(jmax) +
                       kOneTwentyFourth * Power2(jmin)) +
       Power2(amax) * Power2(jmin) *
           (-0.25 * Power2(GetA(s)) + 0.5 * jmax * GetV(s) -
            0.5 * jmax * vtrgt) +
       amax * Power2(jmin) *
           (kOneThird * Power3(GetA(s)) + GetA(s) * jmax * (-GetV(s) + vtrgt) +
            Power2(jmax) * (GetP(s) - ptrgt + texe * vtrgt)) +
       Power2(jmin) *
           (-0.125 * Power4(GetA(s)) +
            Power2(GetA(s)) * jmax * (0.5 * GetV(s) - 0.5 * vtrgt) +
            Power2(jmax) * (-0.5 * Power2(GetV(s)) + GetV(s) * vtrgt -
                            0.5 * Power2(vtrgt)))) /
      (amax * Power2(jmax) * Power2(jmin));

  // postrap hold phase delta velocity
  std::array<double, 3> delta_velocity_hold_phase_postrap_poly_coeff;
  delta_velocity_hold_phase_postrap_poly_coeff[2] = 0.5 * (1 / jmax - 1 / jmin);
  delta_velocity_hold_phase_postrap_poly_coeff[1] = 0.0;
  delta_velocity_hold_phase_postrap_poly_coeff[0] =
      (0.5 * (GetA(s) - amax) * (GetA(s) + amax)) / jmax +
      (0.5 * Power2(amax)) / jmin - GetV(s) + vtrgt;

  // zero time

  std::array<double, 3> zero_time_poly_coeff;
  zero_time_poly_coeff[2] = (0.5 * (jmax - jmin)) / (amax * jmax * jmin);
  zero_time_poly_coeff[1] = 1 / jmax - 1 / jmin;
  zero_time_poly_coeff[0] =
      (0.5 * (Power2(amax) * (jmax - jmin) +
              2. * amax * jmin * (GetA(s) + jmax * texe) -
              jmin * (Power2(GetA(s)) + 2. * jmax * (-GetV(s) + vtrgt)))) /
      (amax * jmax * jmin);

  // Find solution using bisection method
  double peak_acceleration_negtri = GetRootABKMethod(
      absl::bind_front(&GetPositionErrorForPosTrapZeroNegTri,
                       std::cref(position_error_poly_coeff),
                       std::cref(delta_velocity_hold_phase_postrap_poly_coeff),
                       std::cref(zero_time_poly_coeff), std::cref(s)),
      GetNegTriPeakAccelerationLimitsForPosTrapZeroNegTri(
          desired_execution_time, s));

  auto [position_error_optimal, invalid_solution] =
      GetPositionErrorForPosTrapZeroNegTri(
          position_error_poly_coeff,
          delta_velocity_hold_phase_postrap_poly_coeff, zero_time_poly_coeff, s,
          peak_acceleration_negtri);
  if (invalid_solution ||
      (fabs(position_error_optimal) >
       (fabs(kRelPositionErrorTolerance * (GetP(s) - GetPTrgt(s))) +
        kAbsPositionErrorTolerance))) {
    if ((invalid_solution) || (GetTSync(s) < kMaxAllowedSyncTimeSeconds)) {
      return SetFailure(s);
    }
  }
  peak_acceleration_negtri =
      ClipToRange(peak_acceleration_negtri, GetAMin(s), 0.);

  return GeneratePolynomials(
             VUpToVTrgtPosTrapAndThen(HoldToTSyncAndThen(
                 ADownToAPrimeNegTri(peak_acceleration_negtri))),
             s) |
         FinishBlueProfile(s, Profile::kPosTrapZeroNegTri);
}

SyncDOFsMotionState CalcTrajPosTriZeroNegTrap(const SyncDOFsMotionState& s) {
  double desired_execution_time = GetTSync(s) - GetT(s);
  // compute polynomial coefficients

  double jmax = GetJMax(s);

  double amin = GetAMin(s);
  double jmin = GetJMin(s);

  double vtrgt = GetVTrgt(s);
  double ptrgt = GetPTrgt(s);

  double texe = desired_execution_time;

  // postri peak acceleration
  // quartic polynomial coefficients
  // dp = a*x^4 + b*x^3 + c*x^2 + d*x + e
  std::array<double, 5> position_error_poly_coeff;
  position_error_poly_coeff[4] =
      (0.125 * Power2(jmax) - 0.25 * jmax * jmin + 0.125 * Power2(jmin)) /
      (amin * Power2(jmax) * Power2(jmin));
  position_error_poly_coeff[3] = -kOneThird / Power2(jmax) -
                                 kOneSixth / Power2(jmin) + 0.5 / (jmax * jmin);
  position_error_poly_coeff[2] =
      (Power2(GetA(s)) * (0.25 * jmax - 0.25 * jmin) +
       Power2(amin) * (-0.25 * jmax + 0.25 * jmin) +
       GetA(s) * amin * (-0.5 * jmax + 0.5 * jmin) +
       amin * jmax * (-0.5 * jmax + 0.5 * jmin) * texe +
       jmax * (-0.5 * jmax * GetV(s) + 0.5 * jmin * GetV(s) +
               0.5 * jmax * vtrgt - 0.5 * jmin * vtrgt)) /
      (amin * Power2(jmax) * jmin);
  position_error_poly_coeff[1] = 0.0;
  position_error_poly_coeff[0] =
      (Power4(amin) *
           (kOneTwentyFourth * Power2(jmax) - kOneTwentyFourth * Power2(jmin)) +
       amin * Power2(jmin) *
           (-kOneSixth * Power3(GetA(s)) - 0.5 * Power2(GetA(s)) * jmax * texe +
            Power2(jmax) * (GetP(s) - ptrgt + texe * GetV(s))) +
       Power2(amin) * Power2(jmin) *
           (-0.25 * Power2(GetA(s)) + 0.5 * jmax * GetV(s) -
            0.5 * jmax * vtrgt) +
       Power2(jmin) * (0.125 * Power4(GetA(s)) +
                       Power2(GetA(s)) * jmax * (-0.5 * GetV(s) + 0.5 * vtrgt) +
                       Power2(jmax) * (0.5 * Power2(GetV(s)) - GetV(s) * vtrgt +
                                       0.5 * Power2(vtrgt)))) /
      (amin * Power2(jmax) * Power2(jmin));

  // negtrap hold phase delta velocity
  std::array<double, 3> delta_velocity_hold_phase_negtrap_poly_coeff;
  delta_velocity_hold_phase_negtrap_poly_coeff[2] =
      0.5 * (-1. / jmax + 1 / jmin);
  delta_velocity_hold_phase_negtrap_poly_coeff[1] = 0.0;
  delta_velocity_hold_phase_negtrap_poly_coeff[0] =
      (0.5 * (Power2(amin) * (-jmax + jmin) +
              jmin * (Power2(GetA(s)) + 2. * jmax * (-GetV(s) + vtrgt)))) /
      (jmax * jmin);

  // zero time
  std::array<double, 3> zero_time_poly_coeff;
  zero_time_poly_coeff[2] = (0.5 * (-jmax + jmin)) / (amin * jmax * jmin);
  zero_time_poly_coeff[1] = -1. / jmax + 1 / jmin;
  zero_time_poly_coeff[0] =
      (0.5 * (Power2(amin) * (-jmax + jmin) +
              2. * amin * jmin * (GetA(s) + jmax * texe) -
              jmin * (Power2(GetA(s)) + 2. * jmax * (-GetV(s) + vtrgt)))) /
      (amin * jmax * jmin);

  // Find solution using bisection method
  double peak_acceleration_postri = GetRootABKMethod(
      absl::bind_front(&GetPositionErrorForPosTriZeroNegTrap,
                       std::cref(position_error_poly_coeff),
                       std::cref(delta_velocity_hold_phase_negtrap_poly_coeff),
                       std::cref(zero_time_poly_coeff), std::cref(s)),
      GetPosTriPeakAccelerationLimitsForPosTriZeroNegTrap(
          desired_execution_time, s));

  auto [position_error_optimal, invalid_solution] =
      GetPositionErrorForPosTriZeroNegTrap(
          position_error_poly_coeff,
          delta_velocity_hold_phase_negtrap_poly_coeff, zero_time_poly_coeff, s,
          peak_acceleration_postri);
  if (invalid_solution ||
      (fabs(position_error_optimal) >
       (fabs(kRelPositionErrorTolerance * (GetP(s) - GetPTrgt(s))) +
        kAbsPositionErrorTolerance))) {
    if ((invalid_solution) || (GetTSync(s) < kMaxAllowedSyncTimeSeconds)) {
      return SetFailure(s);
    }
  }

  peak_acceleration_postri =
      ClipToRange(peak_acceleration_postri, GetA(s), GetAMax(s));

  return GeneratePolynomials(AUpToAPrimePosTri(peak_acceleration_postri) |
                                 HoldToTSyncAndThen(VDownToVTrgtNegTrap),
                             s) |
         FinishBlueProfile(s, Profile::kPosTriZeroNegTrap);
}

SyncDOFsMotionState CalcTrajPosTriZeroNegTri(const SyncDOFsMotionState& s) {
  double desired_execution_time = GetTSync(s) - GetT(s);

  // Find solution using bisection method
  double peak_acceleration_postri = GetRootABKMethod(
      absl::bind_front(&GetPositionErrorForPosTriZeroNegTri, std::cref(s)),
      GetPosTriPeakAccelerationLimitsForPosTriZeroNegTri(
          s, desired_execution_time));

  auto [position_error_optimal, invalid_solution] =
      GetPositionErrorForPosTriZeroNegTri(s, peak_acceleration_postri);
  if (invalid_solution ||
      (fabs(position_error_optimal) >
       (fabs(kRelPositionErrorTolerance * (GetP(s) - GetPTrgt(s))) +
        kAbsPositionErrorTolerance))) {
    if ((invalid_solution) || (GetTSync(s) < kMaxAllowedSyncTimeSeconds)) {
      return SetFailure(s);
    }
  }
  peak_acceleration_postri =
      ClipToRange(peak_acceleration_postri, GetA(s), GetAMax(s));

  return GeneratePolynomials(AUpToAPrimePosTri(peak_acceleration_postri) |
                                 HoldToTSyncAndThen(VDownToVTrgtNegTri),
                             s) |
         FinishBlueProfile(s, Profile::kPosTriZeroNegTri);
}

SyncDOFsMotionState CalcTrajPosTrapZeroPosTrap(const SyncDOFsMotionState& s) {
  double desired_execution_time = GetTSync(s) - GetT(s);

  // compute polynomial coefficients

  double amax = GetAMax(s);
  double jmax = GetJMax(s);

  double jmin = GetJMin(s);

  double vtrgt = GetVTrgt(s);
  double ptrgt = GetPTrgt(s);

  double texe = desired_execution_time;

  // postrap hold time
  // linear polynomial coefficients
  // dp = a*x + b
  std::array<double, 2> position_error_poly_coeff;
  position_error_poly_coeff[1] =
      (-0.5 * Power2(GetA(s)) * jmin +
       amax * (0.5 * amax * jmax + GetA(s) * jmin - 0.5 * amax * jmin +
               jmax * jmin * texe)) /
          (jmax * jmin) +
      GetV(s) - vtrgt;
  position_error_poly_coeff[0] =
      (Power4(amax) * (-kOneThird * Power2(jmax) + 0.5 * jmax * jmin -
                       kOneSixth * Power2(jmin)) +
       Power3(amax) * (-0.5 * jmax + 0.5 * jmin) * jmin *
           (GetA(s) + jmax * texe) +
       amax * Power2(jmin) *
           (-kOneSixth * Power3(GetA(s)) - 0.5 * Power2(GetA(s)) * jmax * texe +
            Power2(jmax) * (GetP(s) - ptrgt + texe * GetV(s))) +
       Power2(amax) * Power2(jmin) *
           (-0.25 * Power2(GetA(s)) + 0.5 * jmax * GetV(s) -
            0.5 * jmax * vtrgt) +
       Power2(jmin) * (0.125 * Power4(GetA(s)) +
                       Power2(GetA(s)) * jmax * (-0.5 * GetV(s) + 0.5 * vtrgt) +
                       Power2(jmax) * (0.5 * Power2(GetV(s)) - GetV(s) * vtrgt +
                                       0.5 * Power2(vtrgt)))) /
      (amax * Power2(jmax) * Power2(jmin));

  // postrap2 hold phase delta velocity
  // linear polynomial coefficients
  // dv6 = a*x + b
  std::array<double, 2> delta_velocity_hold_phase_postrap2_poly_coeff;
  delta_velocity_hold_phase_postrap2_poly_coeff[1] = -amax;
  delta_velocity_hold_phase_postrap2_poly_coeff[0] =
      (0.5 * Power2(GetA(s))) / jmax + Power2(amax) * (-1. / jmax + 1 / jmin) -
      GetV(s) + vtrgt;

  // Find solution using bisection method
  double hold_time_postrap1 = GetRootABKMethod(
      absl::bind_front(&GetPositionErrorForPosTrapZeroPosTrap,
                       std::cref(position_error_poly_coeff),
                       std::cref(delta_velocity_hold_phase_postrap2_poly_coeff),
                       std::cref(s)),
      GetPosTrapHoldTimeLimitsForPosTrapZeroPosTrap(desired_execution_time, s));

  auto [position_error_optimal, invalid_solution] =
      GetPositionErrorForPosTrapZeroPosTrap(
          position_error_poly_coeff,
          delta_velocity_hold_phase_postrap2_poly_coeff, s, hold_time_postrap1);
  if (invalid_solution ||
      (fabs(position_error_optimal) >
       (fabs(kRelPositionErrorTolerance * (GetP(s) - GetPTrgt(s))) +
        kAbsPositionErrorTolerance))) {
    if ((invalid_solution) || (GetTSync(s) < kMaxAllowedSyncTimeSeconds)) {
      return SetFailure(s);
    }
  }

  hold_time_postrap1 = std::max(hold_time_postrap1, 0.);

  return GeneratePolynomials(AUpToAMax | HoldForDeltaT(hold_time_postrap1) |
                                 AMaxDownToZero |
                                 HoldToTSyncAndThen(VUpToVTrgtPosTrap),
                             s) |
         FinishBlueProfile(s, Profile::kPosTrapZeroPosTrap);
}

SyncDOFsMotionState CalcTrajPosTrapZeroPosTri(const SyncDOFsMotionState& s) {
  double desired_execution_time = GetTSync(s) - GetT(s);
  // compute polynomial coefficients

  double amax = GetAMax(s);
  double jmax = GetJMax(s);

  double jmin = GetJMin(s);

  double vtrgt = GetVTrgt(s);
  double ptrgt = GetPTrgt(s);

  double texe = desired_execution_time;

  // postri peak acceleration
  // quartic polynomial coefficients
  // dp = a*x^4 + b*x^3 + c*x^2 + d*x + e
  std::array<double, 5> position_error_poly_coeff;
  position_error_poly_coeff[4] =
      (-0.125 * Power2(jmax) + 0.25 * jmax * jmin - 0.125 * Power2(jmin)) /
      (amax * Power2(jmax) * Power2(jmin));
  position_error_poly_coeff[3] =
      kOneSixth / Power2(jmax) + kOneThird / Power2(jmin) - 0.5 / (jmax * jmin);
  position_error_poly_coeff[2] =
      (GetA(s) * amax * (0.5 * jmax - 0.5 * jmin) +
       Power2(GetA(s)) * (-0.25 * jmax + 0.25 * jmin) +
       Power2(amax) * (-0.25 * jmax + 0.25 * jmin) +
       amax * jmax * (0.5 * jmax - 0.5 * jmin) * texe +
       jmax * (0.5 * jmax * GetV(s) - 0.5 * jmin * GetV(s) -
               0.5 * jmax * vtrgt + 0.5 * jmin * vtrgt)) /
      (amax * Power2(jmax) * jmin);
  position_error_poly_coeff[1] = 0.0;
  position_error_poly_coeff[0] =
      (Power4(amax) * (-kOneTwentyFourth * Power2(jmax) +
                       kOneTwentyFourth * Power2(jmin)) +
       Power2(amax) * Power2(jmin) *
           (-0.25 * Power2(GetA(s)) + 0.5 * jmax * GetV(s) -
            0.5 * jmax * vtrgt) +
       amax * Power2(jmin) *
           (kOneThird * Power3(GetA(s)) + GetA(s) * jmax * (-GetV(s) + vtrgt) +
            Power2(jmax) * (GetP(s) - ptrgt + texe * vtrgt)) +
       Power2(jmin) *
           (-0.125 * Power4(GetA(s)) +
            Power2(GetA(s)) * jmax * (0.5 * GetV(s) - 0.5 * vtrgt) +
            Power2(jmax) * (-0.5 * Power2(GetV(s)) + GetV(s) * vtrgt -
                            0.5 * Power2(vtrgt)))) /
      (amax * Power2(jmax) * Power2(jmin));

  // postrap hold phase delta velocity
  std::array<double, 3> delta_velocity_hold_phase_postrap_poly_coeff;
  delta_velocity_hold_phase_postrap_poly_coeff[2] =
      0.5 * (-1. / jmax + 1 / jmin);
  delta_velocity_hold_phase_postrap_poly_coeff[1] = 0.0;
  delta_velocity_hold_phase_postrap_poly_coeff[0] =
      (0.5 * (GetA(s) - amax) * (GetA(s) + amax)) / jmax +
      (0.5 * Power2(amax)) / jmin - GetV(s) + vtrgt;

  // zero time

  std::array<double, 3> zero_time_poly_coeff;
  zero_time_poly_coeff[2] = (0.5 * (-jmax + jmin)) / (amax * jmax * jmin);
  zero_time_poly_coeff[1] = -1. / jmax + 1 / jmin;
  zero_time_poly_coeff[0] =
      (0.5 * (Power2(amax) * (jmax - jmin) +
              2. * amax * jmin * (GetA(s) + jmax * texe) -
              jmin * (Power2(GetA(s)) + 2. * jmax * (-GetV(s) + vtrgt)))) /
      (amax * jmax * jmin);

  // Find solution using bisection method
  double peak_acceleration_postri = GetRootABKMethod(
      absl::bind_front(&GetPositionErrorForPosTrapZeroPosTri,
                       std::cref(position_error_poly_coeff),
                       std::cref(delta_velocity_hold_phase_postrap_poly_coeff),
                       std::cref(zero_time_poly_coeff), std::cref(s)),
      GetPosTriPeakAccelerationLimitsForPosTrapZeroPosTri(
          desired_execution_time, s));

  auto [position_error_optimal, invalid_solution] =
      GetPositionErrorForPosTrapZeroPosTri(
          position_error_poly_coeff,
          delta_velocity_hold_phase_postrap_poly_coeff, zero_time_poly_coeff, s,
          peak_acceleration_postri);
  if (invalid_solution ||
      (fabs(position_error_optimal) >
       (fabs(kRelPositionErrorTolerance * (GetP(s) - GetPTrgt(s))) +
        kAbsPositionErrorTolerance))) {
    if ((invalid_solution) || (GetTSync(s) < kMaxAllowedSyncTimeSeconds)) {
      return SetFailure(s);
    }
  }
  peak_acceleration_postri =
      ClipToRange(peak_acceleration_postri, 0., GetAMax(s));

  return GeneratePolynomials(VUpToVTrgtPosTrapAndThen(HoldToTSyncAndThen(
                                 AUpToAPrimePosTri(peak_acceleration_postri))),
                             s) |
         FinishBlueProfile(s, Profile::kPosTrapZeroPosTri);
}

SyncDOFsMotionState CalcTrajPosTriZeroPosTrap(const SyncDOFsMotionState& s) {
  double desired_execution_time = GetTSync(s) - GetT(s);

  // compute polynomial coefficients
  double amax = GetAMax(s);
  double jmax = GetJMax(s);
  double jmin = GetJMin(s);
  double vtrgt = GetVTrgt(s);
  double ptrgt = GetPTrgt(s);
  double texe = desired_execution_time;

  // postri peak acceleration
  // quartic polynomial coefficients
  // dp = a*x^4 + b*x^3 + c*x^2 + d*x + e
  std::array<double, 5> position_error_poly_coeff;
  position_error_poly_coeff[4] =
      (0.125 * Power2(jmax) - 0.25 * jmax * jmin + 0.125 * Power2(jmin)) /
      (amax * Power2(jmax) * Power2(jmin));
  position_error_poly_coeff[3] = -kOneThird / Power2(jmax) -
                                 kOneSixth / Power2(jmin) + 0.5 / (jmax * jmin);
  position_error_poly_coeff[2] =
      (Power2(amax) * jmax * (-0.25 * jmax + 0.25 * jmin) +
       amax * (-0.5 * jmax + 0.5 * jmin) * jmin * (GetA(s) + jmax * texe) +
       jmin * (Power2(GetA(s)) * (0.25 * jmax - 0.25 * jmin) +
               jmax * (-0.5 * jmax * GetV(s) + 0.5 * jmin * GetV(s) +
                       0.5 * jmax * vtrgt - 0.5 * jmin * vtrgt))) /
      (amax * Power2(jmax) * Power2(jmin));
  position_error_poly_coeff[1] = 0.0;
  position_error_poly_coeff[0] =
      (Power4(amax) * (-kOneTwentyFourth * Power2(jmax) +
                       kOneTwentyFourth * Power2(jmin)) +
       amax * Power2(jmin) *
           (-kOneSixth * Power3(GetA(s)) - 0.5 * Power2(GetA(s)) * jmax * texe +
            Power2(jmax) * (GetP(s) - ptrgt + texe * GetV(s))) +
       Power2(amax) * jmax * jmin *
           (-0.25 * Power2(GetA(s)) + 0.5 * jmax * GetV(s) -
            0.5 * jmax * vtrgt) +
       Power2(jmin) * (0.125 * Power4(GetA(s)) +
                       Power2(GetA(s)) * jmax * (-0.5 * GetV(s) + 0.5 * vtrgt) +
                       Power2(jmax) * (0.5 * Power2(GetV(s)) - GetV(s) * vtrgt +
                                       0.5 * Power2(vtrgt)))) /
      (amax * Power2(jmax) * Power2(jmin));

  // postrap hold phase delta velocity
  std::array<double, 3> delta_velocity_hold_phase_postrap_poly_coeff;
  delta_velocity_hold_phase_postrap_poly_coeff[2] =
      0.5 * (-1. / jmax + 1 / jmin);
  delta_velocity_hold_phase_postrap_poly_coeff[1] = 0.0;
  delta_velocity_hold_phase_postrap_poly_coeff[0] =
      (0.5 * (Power2(amax) * (jmax - 1. * jmin) +
              jmin * (Power2(GetA(s)) + 2. * jmax * (-1. * GetV(s) + vtrgt)))) /
      (jmax * jmin);

  // zero time

  std::array<double, 3> zero_time_poly_coeff;
  zero_time_poly_coeff[2] = (0.5 * (-1. * jmax + jmin)) / (amax * jmax * jmin);
  zero_time_poly_coeff[1] = -1. / jmax + 1 / jmin;
  zero_time_poly_coeff[0] =
      (0.5 *
       (Power2(amax) * (jmax - 1. * jmin) +
        2. * amax * jmin * (GetA(s) + jmax * texe) -
        1. * jmin * (Power2(GetA(s)) + 2. * jmax * (-1. * GetV(s) + vtrgt)))) /
      (amax * jmax * jmin);

  // Find solution using bisection method
  double peak_acceleration_postri = GetRootABKMethod(
      absl::bind_front(&GetPositionErrorForPosTriZeroPosTrap,
                       std::cref(position_error_poly_coeff),
                       std::cref(delta_velocity_hold_phase_postrap_poly_coeff),
                       std::cref(zero_time_poly_coeff), std::cref(s)),
      GetPosTriPeakAccelerationLimitsForPosTriZeroPosTrap(
          desired_execution_time, s));

  auto [position_error_optimal, invalid_solution] =
      GetPositionErrorForPosTriZeroPosTrap(
          position_error_poly_coeff,
          delta_velocity_hold_phase_postrap_poly_coeff, zero_time_poly_coeff, s,
          peak_acceleration_postri);
  if (invalid_solution ||
      (fabs(position_error_optimal) >
       (fabs(kRelPositionErrorTolerance * (GetP(s) - GetPTrgt(s))) +
        kAbsPositionErrorTolerance))) {
    if ((invalid_solution) || (GetTSync(s) < kMaxAllowedSyncTimeSeconds)) {
      return SetFailure(s);
    }
  }
  peak_acceleration_postri =
      ClipToRange(peak_acceleration_postri, GetA(s), GetAMax(s));

  return GeneratePolynomials(AUpToAPrimePosTri(peak_acceleration_postri) |
                                 HoldToTSyncAndThen(VUpToVTrgtPosTrap),
                             s) |
         FinishBlueProfile(s, Profile::kPosTriZeroPosTrap);
}

SyncDOFsMotionState CalcTrajPosTriZeroPosTri(const SyncDOFsMotionState& s) {
  double desired_execution_time = GetTSync(s) - GetT(s);

  // Find solution using bisection method
  double peak_acceleration_postri1 =
      GetRootABKMethod(absl::bind_front(&GetPositionErrorForPosTriZeroPosTri,
                                        desired_execution_time, std::cref(s)),
                       GetPosTriPeakAccelerationLimitsForPosTriZeroPosTri(
                           desired_execution_time, s));

  auto [position_error_optimal, invalid_solution] =
      GetPositionErrorForPosTriZeroPosTri(desired_execution_time, s,
                                          peak_acceleration_postri1);
  if (invalid_solution ||
      (fabs(position_error_optimal) >
       (fabs(kRelPositionErrorTolerance * (GetP(s) - GetPTrgt(s))) +
        kAbsPositionErrorTolerance))) {
    if ((invalid_solution) || (GetTSync(s) < kMaxAllowedSyncTimeSeconds)) {
      return SetFailure(s);
    }
  }
  peak_acceleration_postri1 =
      ClipToRange(peak_acceleration_postri1, GetA(s), GetAMax(s));

  return GeneratePolynomials(AUpToAPrimePosTri(peak_acceleration_postri1) |
                                 HoldToTSyncAndThen(VUpToVTrgtPosTri),
                             s) |
         FinishBlueProfile(s, Profile::kPosTriZeroPosTri);
}

SyncDOFsMotionState CalcTrajPosLinHldPosTrap(const SyncDOFsMotionState& s) {
  double desired_execution_time = GetTSync(s) - GetT(s);

  // Find solution using bisection method
  double hold_acceleration = GetRootABKMethod(
      absl::bind_front(&GetPositionErrorForPosLinHldPosTrap,
                       desired_execution_time, std::cref(s)),
      GetHoldAccelerationLimitsForPosLinHldPosTrap(desired_execution_time, s));

  auto [position_error_optimal, invalid_solution] =
      GetPositionErrorForPosLinHldPosTrap(desired_execution_time, s,
                                          hold_acceleration);
  if (invalid_solution ||
      (fabs(position_error_optimal) >
       (fabs(kRelPositionErrorTolerance * (GetP(s) - GetPTrgt(s))) +
        kAbsPositionErrorTolerance))) {
    if ((invalid_solution) || (GetTSync(s) < kMaxAllowedSyncTimeSeconds)) {
      return SetFailure(s);
    }
  }
  hold_acceleration = ClipToRange(hold_acceleration, GetA(s), GetAMax(s));

  // Hld
  double hold_time_ahld = 0.0;
  if (hold_acceleration < GetAMax(s)) {
    hold_time_ahld =
        (GetAMax(s) * desired_execution_time +
         GetA(s) * GetAMax(s) / GetJMax(s) - (GetVTrgt(s) - GetV(s)) -
         0.5 * Power2(GetAMax(s)) * GetInvJMaxSubInvJMin(s) -
         0.5 * Power2(GetA(s)) / GetJMax(s)) /
        (GetAMax(s) - hold_acceleration);
  }
  // check for numerical issues when |max_acc - hold_acc| is a very small value
  if (hold_time_ahld < 0.0 || hold_time_ahld > GetTSync(s)) {
    hold_time_ahld = 0.0;
  }

  return GeneratePolynomials(AUpToAPrime(hold_acceleration) |
                                 HoldForDeltaT(hold_time_ahld) |
                                 VUpToVTrgtPosTrap,
                             s) |
         FinishBlueProfile(s, Profile::kPosLinHldPosTrap);
}

SyncDOFsMotionState CalcTrajPosLinHldPosTri(const SyncDOFsMotionState& s) {
  double desired_execution_time = GetTSync(s) - GetT(s);

  // Find solution using bisection method
  double peak_acceleration_postri =
      GetRootABKMethod(absl::bind_front(&GetPositionErrorForPosLinHldPosTri,
                                        desired_execution_time, std::cref(s)),
                       GetPosTriPeakAccelerationLimitsForPosLinHldPosTri(
                           desired_execution_time, s));

  auto [position_error_optimal, invalid_solution] =
      GetPositionErrorForPosLinHldPosTri(desired_execution_time, s,
                                         peak_acceleration_postri);
  if (invalid_solution ||
      (fabs(position_error_optimal) >
       (fabs(kRelPositionErrorTolerance * (GetP(s) - GetPTrgt(s))) +
        kAbsPositionErrorTolerance))) {
    if ((invalid_solution) || (GetTSync(s) < kMaxAllowedSyncTimeSeconds)) {
      return SetFailure(s);
    }
  }
  peak_acceleration_postri =
      ClipToRange(peak_acceleration_postri, GetA(s), GetAMax(s));

  // hold time
  double hold_time_ahld = desired_execution_time -
                          peak_acceleration_postri * GetInvJMaxSubInvJMin(s) +
                          GetA(s) / GetJMax(s);
  hold_time_ahld = std::max(hold_time_ahld, 0.0);

  // hold acceleration
  double hold_acceleration = GetA(s);
  if (hold_time_ahld > 0) {
    hold_acceleration =
        (GetVTrgt(s) - GetV(s) -
         0.5 * GetInvJMaxSubInvJMin(s) * Power2(peak_acceleration_postri) +
         0.5 * Power2(GetA(s)) / GetJMax(s)) /
        hold_time_ahld;
  }
  hold_acceleration = ClipToRange(hold_acceleration, GetA(s), GetAMax(s));

  return GeneratePolynomials(AUpToAPrime(hold_acceleration) |
                                 HoldForDeltaT(hold_time_ahld) |
                                 AUpToAPrimePosTri(peak_acceleration_postri),
                             s) |
         FinishBlueProfile(s, Profile::kPosLinHldPosTri);
}

SyncDOFsMotionState CalcTrajNegLinHldPosTrap(const SyncDOFsMotionState& s) {
  double desired_execution_time = GetTSync(s) - GetT(s);

  // Find solution using bisection method
  double hold_acceleration = GetRootABKMethod(
      absl::bind_front(&GetPositionErrorForNegLinHldPosTrap,
                       desired_execution_time, std::cref(s)),
      GetHoldAccelerationLimitsForNegLinHldPosTrap(desired_execution_time, s));

  auto [position_error_optimal, invalid_solution] =
      GetPositionErrorForNegLinHldPosTrap(desired_execution_time, s,
                                          hold_acceleration);
  if (invalid_solution ||
      (fabs(position_error_optimal) >
       (fabs(kRelPositionErrorTolerance * (GetP(s) - GetPTrgt(s))) +
        kAbsPositionErrorTolerance))) {
    if ((invalid_solution) || (GetTSync(s) < kMaxAllowedSyncTimeSeconds)) {
      return SetFailure(s);
    }
  }
  hold_acceleration = ClipToRange(hold_acceleration, 0, GetA(s));

  // Hld
  double hold_time_ahld = 0.0;
  if (hold_acceleration < GetAMax(s)) {
    hold_time_ahld = (GetAMax(s) * (desired_execution_time -
                                    GetInvJMaxSubInvJMin(s) *
                                        (GetAMax(s) - hold_acceleration) -
                                    GetA(s) / (-GetJMin(s))) -
                      (GetVTrgt(s) - GetV(s)) +
                      0.5 * GetInvJMaxSubInvJMin(s) *
                          (Power2(GetAMax(s)) - Power2(hold_acceleration)) +
                      0.5 * Power2(GetA(s)) / (-GetJMin(s))) /
                     (GetAMax(s) - hold_acceleration);
  }
  // check for numerical issues when |max_acc - hold_acc| is a very small value
  if (hold_time_ahld < 0.0 || hold_time_ahld > GetTSync(s)) {
    hold_time_ahld = 0.0;
  }

  return GeneratePolynomials(ADownToAPrime(hold_acceleration) |
                                 HoldForDeltaT(hold_time_ahld) |
                                 VUpToVTrgtPosTrap,
                             s) |
         FinishBlueProfile(s, Profile::kNegLinHldPosTrap);
}

SyncDOFsMotionState CalcTrajNegLinHldPosTri(const SyncDOFsMotionState& s) {
  double desired_execution_time = GetTSync(s) - GetT(s);

  // Find solution using bisection method
  double hold_acceleration = GetRootABKMethod(
      absl::bind_front(&GetPositionErrorForNegLinHldPosTri,
                       desired_execution_time, std::cref(s)),
      GetHoldAccelerationLimitsForNegLinHldPosTri(desired_execution_time, s));

  auto [position_error_optimal, invalid_solution] =
      GetPositionErrorForNegLinHldPosTri(desired_execution_time, s,
                                         hold_acceleration);
  if (invalid_solution ||
      (fabs(position_error_optimal) >
       (fabs(kRelPositionErrorTolerance * (GetP(s) - GetPTrgt(s))) +
        kAbsPositionErrorTolerance))) {
    if ((invalid_solution) || (GetTSync(s) < kMaxAllowedSyncTimeSeconds)) {
      return SetFailure(s);
    }
  }
  hold_acceleration = ClipToRange(hold_acceleration, 0, GetA(s));

  // peak acceleration

  std::array<double, 3> coefficients;
  coefficients[2] = 0.5 * GetInvJMaxSubInvJMin(s);
  coefficients[1] = -hold_acceleration * GetInvJMaxSubInvJMin(s);
  coefficients[0] =
      hold_acceleration * (desired_execution_time +
                           0.5 * GetInvJMaxSubInvJMin(s) * hold_acceleration -
                           GetA(s) / (-GetJMin(s))) +
      0.5 * Power2(GetA(s)) / (-GetJMin(s)) - (GetVTrgt(s) - GetV(s));

  PolynomialRoots roots = SafeCalculatePolynomialRoots(coefficients);
  double peak_acceleration_postri = roots[0];
  peak_acceleration_postri =
      ClipToRange(peak_acceleration_postri, hold_acceleration, GetAMax(s));

  // Hld
  double hold_time_ahld =
      desired_execution_time -
      GetInvJMaxSubInvJMin(s) * (peak_acceleration_postri - hold_acceleration) -
      GetA(s) / (-GetJMin(s));
  hold_time_ahld = std::max(hold_time_ahld, 0.0);

  return GeneratePolynomials(ADownToAPrime(hold_acceleration) |
                                 HoldForDeltaT(hold_time_ahld) |
                                 AUpToAPrimePosTri(peak_acceleration_postri),
                             s) |
         FinishBlueProfile(s, Profile::kNegLinHldPosTri);
}

SyncDOFsMotionState CalcTrajPosTrapHldNegLin(const SyncDOFsMotionState& s) {
  double desired_execution_time = GetTSync(s) - GetT(s);

  // Find solution using bisection method
  double hold_time_postrap = GetRootABKMethod(
      absl::bind_front(&GetPositionErrorForPosTrapHldNegLin,
                       desired_execution_time, std::cref(s)),
      GetPosTrapHoldTimeLimitsForPosTrapHldNegLin(desired_execution_time, s));

  auto [position_error_optimal, invalid_solution] =
      GetPositionErrorForPosTrapHldNegLin(desired_execution_time, s,
                                          hold_time_postrap);
  if (invalid_solution ||
      (fabs(position_error_optimal) >
       (fabs(kRelPositionErrorTolerance * (GetP(s) - GetPTrgt(s))) +
        kAbsPositionErrorTolerance))) {
    if ((invalid_solution) || (GetTSync(s) < kMaxAllowedSyncTimeSeconds)) {
      return SetFailure(s);
    }
  }
  hold_time_postrap = std::max(hold_time_postrap, 0.);

  return GeneratePolynomials(
             AUpToAMax | HoldForDeltaT(hold_time_postrap) |
                 ADownToAHldThenHoldToVTrgtTsyncAndThen(ADownToZero),
             s) |
         FinishBlueProfile(s, Profile::kPosTrapHldNegLin);
}

SyncDOFsMotionState CalcTrajPosTriHldNegLin(const SyncDOFsMotionState& s) {
  double desired_execution_time = GetTSync(s) - GetT(s);

  // Find solution using bisection method
  double peak_acceleration_postri =
      GetRootABKMethod(absl::bind_front(&GetPositionErrorForPosTriHldNegLin,
                                        desired_execution_time, std::cref(s)),
                       GetPosTriPeakAccelerationLimitsForPosTriHldNegLin(
                           desired_execution_time, s));

  auto [position_error_optimal, invalid_solution] =
      GetPositionErrorForPosTriHldNegLin(desired_execution_time, s,
                                         peak_acceleration_postri);
  if (invalid_solution ||
      (fabs(position_error_optimal) >
       (fabs(kRelPositionErrorTolerance * (GetP(s) - GetPTrgt(s))) +
        kAbsPositionErrorTolerance))) {
    if ((invalid_solution) || (GetTSync(s) < kMaxAllowedSyncTimeSeconds)) {
      return SetFailure(s);
    }
  }
  peak_acceleration_postri =
      ClipToRange(peak_acceleration_postri, GetA(s), GetAMax(s));

  return GeneratePolynomials(
             AUpToAPrime(peak_acceleration_postri) |
                 ADownToAHldThenHoldToVTrgtTsyncAndThen(ADownToZero),
             s) |
         FinishBlueProfile(s, Profile::kPosTriHldNegLin);
}

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic
