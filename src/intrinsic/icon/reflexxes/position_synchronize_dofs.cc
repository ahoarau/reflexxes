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

// This file defines the decision tree for the Position OTG algorithm
// SynchronizeDOFs function otherwise known as "Position Step2" in the
// Reflexxes literature.
// See
// https://github.com/intrinsic-ai/intrinsic-core/blob/main/intrinsic_control/intrinsic/icon/reflexxes/g3doc/step2_decision_tree.pdf
// for a detailed description of what's happening in this file.

#include <algorithm>
#include <array>

#include "intrinsic/icon/reflexxes/inputs.h"
#include "intrinsic/icon/reflexxes/internal/decision_tree_utility_functions.h"
#include "intrinsic/icon/reflexxes/internal/intermediate_profile_funcs.h"
#include "intrinsic/icon/reflexxes/internal/math.h"
#include "intrinsic/icon/reflexxes/internal/motion_state.h"
#include "intrinsic/icon/reflexxes/internal/motion_state_transforms.h"
#include "intrinsic/icon/reflexxes/internal/polynomial_solvers.h"
#include "intrinsic/icon/reflexxes/internal/position_calc_traj_funcs.h"
#include "intrinsic/icon/reflexxes/internal/synchronize_dofs_base.h"
#include "intrinsic/icon/reflexxes/position_flags.h"
#include "intrinsic/icon/reflexxes/position_outputs.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {
namespace {

// if a->amax->hold->0 and then a->+apeak->0, so that v=vtrgt and t=tsync, is
// p<=ptrgt?
bool CalcDecision16(const SyncDOFsMotionState& s) {
  // time to make a->amax using jmax
  double acceleration_time_amax = (GetAMax(s) - GetA(s)) / GetJMax(s);

  // change in velocity during acceleration phase to amax
  double delta_velocity_acceleration_phase_amax =
      0.5 * (Power2(GetAMax(s)) - Power2(GetA(s))) / GetJMax(s);

  // time to make a(=amax)->0 using jmin
  double deceleration_time_amax = GetAMax(s) / (-GetJMin(s));

  // change in velocity during deceleration phase from amax
  double delta_velocity_deceleration_phase_amax =
      0.5 * Power2(GetAMax(s)) / (-GetJMin(s));

  // desired change in velocity
  double delta_velocity_desired = GetVTrgt(s) - GetV(s);

  // desired execution time
  double delta_time_desired = GetTSync(s) - GetT(s);

  // constant value (this needs to be negative for a valid solution
  double const_value =
      (2 / GetInvJMaxSubInvJMin(s)) *
      (delta_velocity_desired - delta_velocity_acceleration_phase_amax -
       delta_velocity_deceleration_phase_amax -
       GetAMax(s) * (delta_time_desired - acceleration_time_amax -
                     deceleration_time_amax));

  // peak acceleration
  double peak_acceleration =
      GetAMax(s) - GetSqrt(Power2(GetAMax(s)) + const_value);

  return AUpToAMax(s) |
         HoldToTSyncAndThen(AMaxDownToZero |
                            AUpToAPrimePosTri(peak_acceleration)) |
         CmpPLtePTrgt;
}

// if a-up->+ahld->hold->0 so that t=tsync and v=vtrgt, is p<=ptrgt?
bool CalcDecision22(const SyncDOFsMotionState& s) {
  // Figure out the hold acceleration.
  std::array<double, 3> coefficients;
  coefficients[2] = 0.5 * GetInvJMaxSubInvJMin(s);
  coefficients[1] = -(GetTSync(s) - GetT(s)) - GetA(s) / GetJMax(s);
  coefficients[0] =
      0.5 * Power2(GetA(s)) / GetJMax(s) + (GetVTrgt(s) - GetV(s));

  PolynomialRoots roots = SafeCalculatePolynomialRoots(coefficients);

  double hold_acceleration = GetA(s);
  for (double root : roots) {
    if (root > GetA(s) && root <= GetAMax(s) &&
        (GetTSync(s) - GetT(s) - root * GetInvJMaxSubInvJMin(s) +
         GetA(s) / GetJMax(s)) >= 0) {
      hold_acceleration = root;
    }
  }

  return AUpToAPrime(hold_acceleration, s) | HoldToTSyncAndThen(ADownToZero) |
         CmpPLtePTrgt;
}

// if a->hold and then v->vtrgt(PosTri) so that t=tsync, is p<=ptrgt?
bool CalcDecision24(const SyncDOFsMotionState& s) {
  // Figure out peak acceleration.
  std::array<double, 3> coefficients;
  coefficients[2] = 0.5 * GetInvJMaxSubInvJMin(s);
  coefficients[1] = -GetA(s) * GetInvJMaxSubInvJMin(s);
  coefficients[0] = 0.5 * Power2(GetA(s)) / GetJMax(s) -
                    (GetVTrgt(s) - GetV(s)) + GetA(s) * (GetTSync(s) - GetT(s));
  double peak_acceleration = SafeCalculatePolynomialRoots(coefficients)[0];
  peak_acceleration = ClipToRange(peak_acceleration, GetA(s), GetAMax(s));

  return HoldToTSyncAndThen(AUpToAPrimePosTri(peak_acceleration), s) |
         CmpPLtePTrgt;
}

// if a->hold and then v->vtrgt(PosTrap) so that t=tsync, is p<=ptrgt?
bool CalcDecision25(const SyncDOFsMotionState& s) {
  // hold time at a
  double hold_time_a =
      (GetAMax(s) * (GetTSync(s) - GetT(s)) - (GetVTrgt(s) - GetV(s)) -
       0.5 * Power2(GetAMax(s)) * GetInvJMaxSubInvJMin(s) +
       GetA(s) * GetAMax(s) / GetJMax(s) - 0.5 * Power2(GetA(s)) / GetJMax(s)) /
      (GetAMax(s) - GetA(s));

  return HoldForDeltaT(hold_time_a, s) | VUpToVTrgtPosTrap | CmpPLtePTrgt;
}

// if a->+apeak->0 and then v->vtrgt(PosTrap) so that t=tsync, is p<=ptrgt?
bool CalcDecision31(const SyncDOFsMotionState& s) {
  // Calc peak acceleration.
  std::array<double, 3> coefficients;
  coefficients[2] = 0.5 * GetInvJMaxSubInvJMin(s);
  coefficients[1] = -GetAMax(s) * GetInvJMaxSubInvJMin(s);
  coefficients[0] =
      GetAMax(s) * (GetTSync(s) - GetT(s)) + GetA(s) * GetAMax(s) / GetJMax(s) -
      0.5 * Power2(GetAMax(s)) * GetInvJMaxSubInvJMin(s) -
      0.5 * Power2(GetA(s)) / GetJMax(s) - (GetVTrgt(s) - GetV(s));

  double peak_acceleration = SafeCalculatePolynomialRoots(coefficients)[1];
  peak_acceleration = std::max(peak_acceleration, GetA(s));

  return AUpToAPrimePosTri(peak_acceleration, s) | VUpToVTrgtPosTrap |
         CmpPLtePTrgt;
}

// if a->+apeak1->0 and then a->+apeak2->0 so that v=vtrgt (either solution), is
// t>=tsync?
bool CalcDecision37(const SyncDOFsMotionState& s) {
  // peak acceleration if apeak1 = apeak2 = apeak
  double peak_acceleration = GetSqrt(
      0.5 * (2.0 * (GetVTrgt(s) - GetV(s)) + Power2(GetA(s)) / GetJMax(s)) /
      GetInvJMaxSubInvJMin(s));

  // time durations of the acceleration and deceleration phases of a->apeak1->0
  double acceleration_time_apeak1 = (peak_acceleration - GetA(s)) / GetJMax(s);
  double deceleration_time_apeak1 = peak_acceleration / (-GetJMin(s));

  // time durations of the acceleration and deceleration phases of a->apeak2->0
  double acceleration_time_apeak2 = peak_acceleration / GetJMax(s);
  double deceleration_time_apeak2 = peak_acceleration / (-GetJMin(s));

  if ((GetT(s) + acceleration_time_apeak1 + deceleration_time_apeak1 +
       acceleration_time_apeak2 + deceleration_time_apeak2) < GetTSync(s)) {
    return false;
  }

  // peak acceleration1 - max and min values
  std::array<double, 3> coefficients;
  coefficients[2] = 1;
  coefficients[1] = -1.0 * (GetTSync(s) - GetT(s) + (GetA(s) / GetJMax(s))) /
                    GetInvJMaxSubInvJMin(s);
  coefficients[0] =
      (0.5 * Power2((GetTSync(s) - GetT(s) + (GetA(s) / GetJMax(s))) /
                    GetInvJMaxSubInvJMin(s))) -
      (0.5 * (2.0 * (GetVTrgt(s) - GetV(s)) + (Power2(GetA(s)) / GetJMax(s))) /
       GetInvJMaxSubInvJMin(s));
  PolynomialRoots roots = SafeCalculatePolynomialRoots(coefficients);

  double peak1_acceleration_max = roots[0];
  peak1_acceleration_max = std::max(peak1_acceleration_max, GetA(s));
  if (AUpToAPrimePosTri(peak1_acceleration_max, s) | VUpToVTrgtPosTri |
      CmpPGtPTrgt) {
    return true;
  }

  double peak1_acceleration_min = roots[1];
  peak1_acceleration_min = std::max(peak1_acceleration_min, GetA(s));

  if (AUpToAPrimePosTri(peak1_acceleration_min, s) | VUpToVTrgtPosTri |
      CmpPLtPTrgt) {
    return true;
  }

  return false;
}

// The decision tree.
// This is done as a struct to allow the decisions to be listed in order
// starting from the first decision.  Bundling the decisions into a struct
// allows us to refer to later decisions before they are declared, whereas a
// namespace wouldn't allow that.
struct Decisions {
  inline static SyncDOFsMotionState D1(const SyncDOFsMotionState& s) {
    return MakeDecision(1, CmpAGteZero, &D2, FlipState | &D2, s);
  }

  inline static SyncDOFsMotionState D2(const SyncDOFsMotionState& s) {
    return MakeDecision(2, CmpALteAMax, &D3, CalcTrajNegLinAdowntoAMax | &D3,
                        s);
  }

  inline static SyncDOFsMotionState D3(const SyncDOFsMotionState& s) {
    return MakeDecision(3, ADownToZero | CmpVLteVMax, &D4,
                        CalcTrajNegLinAdowntoZero | FlipState |
                            IfElse(GetIntoBoundariesFast(s), &D5, &D72),
                        s);
  }

  inline static SyncDOFsMotionState D4(const SyncDOFsMotionState& s) {
    return MakeDecision(4, CmpVGteVMin, &D9,
                        IfElse(GetIntoBoundariesFast(s), &D5, &D72), s);
  }

  inline static SyncDOFsMotionState D5(const SyncDOFsMotionState& s) {
    return MakeDecision(5, AUpToAMax | CmpVGtVMin, &D6, &D7, s);
  }

  inline static SyncDOFsMotionState D6(const SyncDOFsMotionState& s) {
    return MakeDecision(6, VUpToVMinPosLin | ADownToZero | CmpVLteVMax,
                        CalcTrajPosLinVuptoVMin | &D9,
                        CalcTrajPosLinNegLinVuptoVMin | &D9, s);
  }

  inline static SyncDOFsMotionState D7(const SyncDOFsMotionState& s) {
    // if v->vmin (PosLinHld) and then a->0, is v>vmax?
    // We just assume the end of vmin poslin
    return MakeDecision(
        7, SetV(GetVMin(s)) | SetA(GetAMax(s)) | ADownToZero | CmpVGteVMax, &D8,
        CalcTrajPosLinHldVuptoVMin | &D9, s);
  }

  inline static SyncDOFsMotionState D8(const SyncDOFsMotionState& s) {
    return MakeDecision(8, AUpToAMaxDownToZero | CmpVGtVMax,
                        CalcTrajPosLinNegLinVuptoVMin | &D9,
                        CalcTrajPosLinHldNegLinVuptoVMin | &D9, s);
  }

  inline static SyncDOFsMotionState D9(const SyncDOFsMotionState& s) {
    return MakeDecision(9, ADownToZero | CmpVLteVTrgt, &D10, &D62, s);
  }

  inline static SyncDOFsMotionState D10(const SyncDOFsMotionState& s) {
    return MakeDecision(10, AUpToAMaxDownToZero | CmpVLteVTrgt, &D11, &D50, s);
  }

  inline static SyncDOFsMotionState D11(const SyncDOFsMotionState& s) {
    // if v->vtrgt(PosTrap)->hold so that t=tsync, is p<=ptrgt?
    return MakeDecision(
        11, VUpToVTrgtPosTrap | HoldToTSync | CmpPLtePTrgtVMaxEpsilon, &D12,
        &D14, s);
  }

  inline static SyncDOFsMotionState D12(const SyncDOFsMotionState& s) {
    // if a->amax->hold->0 and then a->amin->0, so that v=vtrgt, is t>=tsync?
    return MakeDecision(
        12,
        AUpToAMax | HoldToVTrgtAndThen(ADownToZero | ADownToAMinUpToZero) |
            CmpTGteTSync,
        CalcTrajPosTrapZeroNegTri, &D13, s);
  }

  inline static SyncDOFsMotionState D13(const SyncDOFsMotionState& s) {
    // if a->amax->hold->0, then a(=0)->hold and then a->amin->0, so that
    // t=tsync
    // and v=vtrgt, is p<=ptrgt?
    return MakeDecision(
        13,
        AUpToAMax |
            HoldToVTrgtAndThen(ADownToZero |
                               HoldToTSyncAndThen(ADownToAMinUpToZero)) |
            CmpPLtePTrgt,
        CalcTrajPosTrapZeroNegTrap, CalcTrajPosTrapZeroNegTri, s);
  }

  inline static SyncDOFsMotionState D14(const SyncDOFsMotionState& s) {
    // if a->0 and then a->amax->0, is v<=vtrgt?
    return MakeDecision(14, ADownToZero | AUpToAMaxDownToZero | CmpVLteVTrgt,
                        &D15, &D41, s);
  }

  inline static SyncDOFsMotionState D15(const SyncDOFsMotionState& s) {
    // if a->0 and then v->vtrgt(PosTrap), is t>=tsync?
    return MakeDecision(15, ADownToZero | VUpToVTrgtPosTrap | CmpTGteTSync,
                        &D16, &D27, s);
  }

  inline static SyncDOFsMotionState D16(const SyncDOFsMotionState& s) {
    return MakeDecision(16, &CalcDecision16, CalcTrajPosTrapZeroPosTri, &D17,
                        s);
  }

  inline static SyncDOFsMotionState D17(const SyncDOFsMotionState& s) {
    // if a->hold->0 so that t=tsync, is v>vtrgt?
    return MakeDecision(17, HoldToTSyncAndThen(ADownToZero) | CmpVGtVTrgt, &D18,
                        &D22, s);
  }

  inline static SyncDOFsMotionState D18(const SyncDOFsMotionState& s) {
    // if a->+ahld->hold->0 so that t=tsync and v=vtrgt, is p(epsilon)<=ptrgt?
    return MakeDecision(18,
                        ADownToAHldThenHoldToVTrgtTsyncAndThen(ADownToZero) |
                            CmpPLteEpsilonPTrgt,
                        &D19, &D20, s);
  }

  inline static SyncDOFsMotionState D19(const SyncDOFsMotionState& s) {
    // if a->amax->+ahld->hold->0 so that t=tsync and v=vtrgt, is p<=ptrgt?
    return MakeDecision(
        19,
        AUpToAMax | ADownToAHldThenHoldToVTrgtTsyncAndThen(ADownToZero) |
            CmpPLtePTrgt,
        CalcTrajPosTrapHldNegLin, CalcTrajPosTriHldNegLin, s);
  }

  inline static SyncDOFsMotionState D20(const SyncDOFsMotionState& s) {
    // if a->hold->amax->0 so that t=tsync, is v>vtrgt?
    return MakeDecision(20,
                        HoldToTSyncAndThen(AUpToAMaxDownToZero) | CmpVGtVTrgt,
                        &D70, CalcTrajNegLinHldPosTrap, s);
  }

  inline static SyncDOFsMotionState D22(const SyncDOFsMotionState& s) {
    return MakeDecision(22, &CalcDecision22, &D19, &D23, s);
  }

  inline static SyncDOFsMotionState D23(const SyncDOFsMotionState& s) {
    // if a->hold->amax->0 so that t=tsync, is v>vtrgt?
    return MakeDecision(23,
                        HoldToTSyncAndThen(AUpToAMaxDownToZero) | CmpVGtVTrgt,
                        &D24, &D25, s);
  }

  inline static SyncDOFsMotionState D24(const SyncDOFsMotionState& s) {
    return MakeDecision(24, &CalcDecision24, CalcTrajPosLinHldPosTri, &D70, s);
  }

  inline static SyncDOFsMotionState D25(const SyncDOFsMotionState& s) {
    return MakeDecision(25, &CalcDecision25, &D71, CalcTrajNegLinHldPosTrap, s);
  }

  inline static SyncDOFsMotionState D27(const SyncDOFsMotionState& s) {
    // if a->0 then a->hold and then v->vtrgt(PosTrap) so that t=tsync, is
    // p<=ptrgt?
    return MakeDecision(
        27, ADownToZero | HoldToTSyncAndThen(VUpToVTrgtPosTrap) | CmpPLtePTrgt,
        &D28, CalcTrajNegLinAdowntoZero | FlipState | &D39, s);
  }

  inline static SyncDOFsMotionState D28(const SyncDOFsMotionState& s) {
    // if a->amax->0 and then a->amax->0, is v<=vtrgt?
    return MakeDecision(
        28, AUpToAMaxDownToZero | AUpToAMaxDownToZero | CmpVLteVTrgt, &D29,
        &D34, s);
  }

  inline static SyncDOFsMotionState D29(const SyncDOFsMotionState& s) {
    // if a->amax->0 and then v->vtrgt(PosTrap), is t>=tsync?
    return MakeDecision(29,
                        AUpToAMaxDownToZero | VUpToVTrgtPosTrap | CmpTGteTSync,
                        &D30, &D32, s);
  }

  inline static SyncDOFsMotionState D30(const SyncDOFsMotionState& s) {
    // if a->amax->hold->0 and then a->+apeak->0, so that v=vtrgt and t=tsync,
    // is p<=ptrgt?
    return MakeDecision(30, &CalcDecision16, CalcTrajPosTrapZeroPosTri, &D31,
                        s);
  }

  inline static SyncDOFsMotionState D31(const SyncDOFsMotionState& s) {
    return MakeDecision(31, &CalcDecision31, &D17, CalcTrajPosTriZeroPosTrap,
                        s);
  }

  inline static SyncDOFsMotionState D32(const SyncDOFsMotionState& s) {
    // if a->amax->0, then a(=0)->hold and then v->vtrgt(PosTrap) so that
    // t=tsync,
    // is p<=ptrgt?
    return MakeDecision(32,
                        AUpToAMaxDownToZero |
                            HoldToTSyncAndThen(VUpToVTrgtPosTrap) |
                            CmpPLtePTrgt,
                        &D33, CalcTrajPosTriZeroPosTrap, s);
  }

  inline static SyncDOFsMotionState D33(const SyncDOFsMotionState& s) {
    // if a->amax->hold->0, then a(=0)->hold and then a(=0)->amax->0 so that
    // v=vtrgt
    // and t=tsync, is p<=ptrgt?
    return MakeDecision(
        33,
        AUpToAMax |
            HoldToVTrgtAndThen(ADownToZero |
                               HoldToTSyncAndThen(AUpToAMaxDownToZero)) |
            CmpPLtePTrgt,
        CalcTrajPosTrapZeroPosTri, CalcTrajPosTrapZeroPosTrap, s);
  }

  inline static SyncDOFsMotionState D34(const SyncDOFsMotionState& s) {
    // if a->amax->0 and then v->vtrgt(PosTri), is t>=tsync?
    return MakeDecision(34,
                        AUpToAMaxDownToZero | VUpToVTrgtPosTri | CmpTGteTSync,
                        &D35, &D38, s);
  }

  inline static SyncDOFsMotionState D35(const SyncDOFsMotionState& s) {
    // if a->apeak->0 and then a->amax->0 so that v=vtrgt, is t>=tsync?
    return MakeDecision(
        35, VToVTrgtPosTriAndThen(AUpToAMaxDownToZero) | CmpTGteTSync, &D16,
        &D36, s);
  }

  inline static SyncDOFsMotionState D36(const SyncDOFsMotionState& s) {
    // if a->+apeak->0, then a(=0)->hold and then a(=0)->amax->0 so that v=vtrgt
    // and
    // t=tsync, is p<=ptrgt?
    return MakeDecision(
        36,
        VToVTrgtPosTriAndThen(HoldToTSyncAndThen(AUpToAMaxDownToZero)) |
            CmpPLtePTrgt,
        &D37, CalcTrajPosTriZeroPosTrap, s);
  }

  inline static SyncDOFsMotionState D37(const SyncDOFsMotionState& s) {
    return MakeDecision(37, &CalcDecision37, &D17, CalcTrajPosTriZeroPosTri, s);
  }

  inline static SyncDOFsMotionState D38(const SyncDOFsMotionState& s) {
    // if a->amax->0, then a(=0)->hold and then v->vtrgt(PosTri) so that
    // t=tsync, is
    // p<=ptrgt?
    return MakeDecision(38,
                        AUpToAMaxDownToZero |
                            HoldToTSyncAndThen(VUpToVTrgtPosTri) | CmpPLtePTrgt,
                        CalcTrajPosTrapZeroPosTri, &D35, s);
  }

  inline static SyncDOFsMotionState D39(const SyncDOFsMotionState& s) {
    // if a->amax->0 and then v->vtrgt(NegTrap), is t>=tsync?
    return MakeDecision(
        39, AUpToAMaxDownToZero | VDownToVTrgtNegTrap | CmpTGteTSync,
        CalcTrajPosTriZeroNegTrap, &D40, s);
  }

  inline static SyncDOFsMotionState D40(const SyncDOFsMotionState& s) {
    // if a->amax->0, then a(=0)->hold and then v->vtrgt(NegTrap) so that
    // t=tsync,
    // is p<=ptrgt?
    return MakeDecision(
        40,
        AUpToAMaxDownToZero | HoldToTSyncAndThen(VDownToVTrgtNegTrap) |
            CmpPLtePTrgt,
        CalcTrajPosTrapZeroNegTrap, CalcTrajPosTriZeroNegTrap, s);
  }

  inline static SyncDOFsMotionState D41(const SyncDOFsMotionState& s) {
    // if a->0 and then v->vtrgt(PosTri), is t>=tsync?
    return MakeDecision(41, ADownToZero | VUpToVTrgtPosTri | CmpTGteTSync, &D42,
                        &D44, s);
  }

  inline static SyncDOFsMotionState D42(const SyncDOFsMotionState& s) {
    // if a->amax->0 and then v->vtrgt(PosTri), is t>=tsync?
    return MakeDecision(42,
                        AUpToAMaxDownToZero | VUpToVTrgtPosTri | CmpTGteTSync,
                        &D16, &D43, s);
  }

  inline static SyncDOFsMotionState D43(const SyncDOFsMotionState& s) {
    // if a->amax->0, then a(=0)->hold and then v->vtrgt(PosTri) so that
    // t=tsync, is
    // p<=ptrgt?
    return MakeDecision(43,
                        AUpToAMaxDownToZero |
                            HoldToTSyncAndThen(VUpToVTrgtPosTri) | CmpPLtePTrgt,
                        CalcTrajPosTrapZeroPosTri, &D17, s);
  }

  inline static SyncDOFsMotionState D44(const SyncDOFsMotionState& s) {
    // if a->0, then a(=0)->hold and then v->vtrgt(PosTri) so that t=tsync, is
    // p(epsilon)<=ptrgt?
    return MakeDecision(44,
                        ADownToZero | HoldToTSyncAndThen(VUpToVTrgtPosTri) |
                            CmpPLteEpsilonPTrgt,
                        &D45, CalcTrajNegLinAdowntoZero | FlipState | &D67, s);
  }

  inline static SyncDOFsMotionState D45(const SyncDOFsMotionState& s) {
    // if a->amax->0 and then v->vtrgt(PosTri), is t>=tsync?
    return MakeDecision(45,
                        AUpToAMaxDownToZero | VUpToVTrgtPosTri | CmpTGteTSync,
                        &D46, &D47, s);
  }

  inline static SyncDOFsMotionState D46(const SyncDOFsMotionState& s) {
    // set ax=a, if a->0 and then a->+apeak->0 so that v=vtrgt, is ax>=apeak?
    return MakeDecision(46,
                        CmpGte(GetA, ADownToZero | CalcPeakAccelVTrgtPosTri),
                        CalcTrajPosTriZeroPosTri, &D37, s);
  }

  inline static SyncDOFsMotionState D47(const SyncDOFsMotionState& s) {
    // if a->amax->0, then a(=0)->hold and then v->vtrgt(PosTri) so that
    // t=tsync, is p<=ptrgt?
    return MakeDecision(47,
                        AUpToAMaxDownToZero |
                            HoldToTSyncAndThen(VUpToVTrgtPosTri) | CmpPLtePTrgt,
                        CalcTrajPosTrapZeroPosTri, &D46, s);
  }

  inline static SyncDOFsMotionState D48(const SyncDOFsMotionState& s) {
    // if a->+apeak->0 and then a->amin->0 so that v=vtrgt, is t>=tsync?
    return MakeDecision(
        48, VToVTrgtPosTriAndThen(ADownToAMinUpToZero) | CmpTGteTSync,
        CalcTrajPosTriZeroNegTri, &D49, s);
  }

  inline static SyncDOFsMotionState D49(const SyncDOFsMotionState& s) {
    // if a->+apeak->0, then a(=0)->hold and then a(=0)->amin->0 so that v=vtrgt
    // and
    // t=tsync, is p<=ptrgt?
    return MakeDecision(
        49,
        VToVTrgtPosTriAndThen(HoldToTSyncAndThen(ADownToAMinUpToZero)) |
            CmpPLtePTrgt,
        &D39, CalcTrajPosTriZeroNegTri, s);
  }

  inline static SyncDOFsMotionState D50(const SyncDOFsMotionState& s) {
    // if v->vtrgt(PosTri)->hold so that t=tsync, is p<=ptrgt?
    return MakeDecision(
        50, VUpToVTrgtPosTri | HoldToTSync | CmpPLtePTrgtVMaxEpsilon, &D67,
        &D53, s);
  }

  inline static SyncDOFsMotionState D51(const SyncDOFsMotionState& s) {
    // if a->amax->0 and then v->vtrgt(NegTri), is t>=tsync?
    return MakeDecision(51,
                        AUpToAMaxDownToZero | VDownToVTrgtNegTri | CmpTGteTSync,
                        CalcTrajPosTriZeroNegTri, &D52, s);
  }

  inline static SyncDOFsMotionState D52(const SyncDOFsMotionState& s) {
    // if a->amax->0, then a(=0)->hold and then v->vtrgt(NegTri) so that
    // t=tsync, is
    // p<=ptrgt?
    return MakeDecision(52,
                        AUpToAMaxDownToZero |
                            HoldToTSyncAndThen(VDownToVTrgtNegTri) |
                            CmpPLtePTrgt,
                        &D12, CalcTrajPosTriZeroNegTri, s);
  }

  inline static SyncDOFsMotionState D53(const SyncDOFsMotionState& s) {
    // if a->0 and then v->vtrgt(PosTri), is t>=tsync?
    return MakeDecision(53, ADownToZero | VUpToVTrgtPosTri | CmpTGteTSync, &D54,
                        &D60, s);
  }

  inline static SyncDOFsMotionState D54(const SyncDOFsMotionState& s) {
    // set ax=a, if a->0 and then a->+apeak->0 so that v=vtrgt, is ax>=apeak?
    return MakeDecision(54,
                        CmpGte(GetA, ADownToZero | CalcPeakAccelVTrgtPosTri),
                        &D55, &D59, s);
  }

  inline static SyncDOFsMotionState D55(const SyncDOFsMotionState& s) {
    // if a->hold->0 so that t=tsync, is v>vtrgt?
    return MakeDecision(55, HoldToTSyncAndThen(ADownToZero) | CmpVGtVTrgt, &D56,
                        &D57, s);
  }

  inline static SyncDOFsMotionState D56(const SyncDOFsMotionState& s) {
    // if a->+ahld->hold->0 so that t=tsync and v=vtrgt, is p<=ptrgt?
    return MakeDecision(56,
                        ADownToAHldThenHoldToVTrgtTsyncAndThen(ADownToZero) |
                            CmpPLteEpsilonPTrgt,
                        CalcTrajPosTriHldNegLin, CalcTrajNegLinHldPosTri, s);
  }

  inline static SyncDOFsMotionState D57(const SyncDOFsMotionState& s) {
    // a-up->+ahld->hold->0 so that t=tsync and v=vtrgt, is p<=ptrgt?
    return MakeDecision(57, &CalcDecision22, CalcTrajPosTriHldNegLin, &D58, s);
  }

  inline static SyncDOFsMotionState D58(const SyncDOFsMotionState& s) {
    // if a->hold and then v->vtrgt(PosTri) so that t=tsync, is p<=ptrgt?
    return MakeDecision(58, &CalcDecision24, CalcTrajPosLinHldPosTri,
                        CalcTrajNegLinHldPosTri, s);
  }

  inline static SyncDOFsMotionState D59(const SyncDOFsMotionState& s) {
    // if a->+apeak1->0 and then a->+apeak2->0 so that v=vtrgt (either
    // solution), is
    // t>=tsync?

    return MakeDecision(59, &CalcDecision37, &D55, CalcTrajPosTriZeroPosTri, s);
  }

  inline static SyncDOFsMotionState D60(const SyncDOFsMotionState& s) {
    // if a->0, then a(=0)->hold and then v->vtrgt(PosTri) so that t=tsync, is
    // p(epsilon)<=ptrgt?
    return MakeDecision(60,
                        ADownToZero | HoldToTSyncAndThen(VUpToVTrgtPosTri) |
                            CmpPLteEpsilonPTrgt,
                        &D61, CalcTrajNegLinAdowntoZero | FlipState | &D67, s);
  }

  inline static SyncDOFsMotionState D61(const SyncDOFsMotionState& s) {
    // set ax=a, if a->0 and then a->+apeak->0 so that v=vtrgt, is ax>=apeak?
    return MakeDecision(61,
                        CmpGte(GetA, ADownToZero | CalcPeakAccelVTrgtPosTri),
                        CalcTrajPosTriZeroPosTri, &D59, s);
  }

  inline static SyncDOFsMotionState D62(const SyncDOFsMotionState& s) {
    // if a->0 and then a->amin->0, is v<=vtrgt?
    return MakeDecision(62, ADownToZero | ADownToAMinUpToZero | CmpVLteVTrgt,
                        &D63, &D68, s);
  }

  inline static SyncDOFsMotionState D63(const SyncDOFsMotionState& s) {
    // if a->0, then v->vtrgt(NegTri) and then a(=0)->hold so that t=tsync, is
    // p>=ptrgt?
    return MakeDecision(63,
                        ADownToZero | VDownToVTrgtNegTri | HoldToTSync |
                            CmpPGtePTrgtVMinEpsilon,
                        CalcTrajNegLinAdowntoZero | FlipState | &D67, &D64, s);
  }

  inline static SyncDOFsMotionState D64(const SyncDOFsMotionState& s) {
    // if a->0, then a(=0)->hold and then v->vtrgt(NegTri) so that t=tsync, is
    // p>=(epsilon)ptrgt?
    return MakeDecision(64,
                        ADownToZero | HoldToTSyncAndThen(VDownToVTrgtNegTri) |
                            CmpPGteEpsilonPTrgt,
                        CalcTrajNegLinAdowntoZero | FlipState | &D65, &D67, s);
  }

  inline static SyncDOFsMotionState D65(const SyncDOFsMotionState& s) {
    // if a->+apeak1->0 and then a->+apeak2->0 so that v=vtrgt (either
    // solution), is
    // t>=tsync?
    return MakeDecision(65, &CalcDecision37, &D66, CalcTrajPosTriZeroPosTri, s);
  }

  inline static SyncDOFsMotionState D66(const SyncDOFsMotionState& s) {
    // if a-up->+ahld->hold->0 so that t=tsync and v=vtrgt, is p<=ptrgt?
    return MakeDecision(66, &CalcDecision22, CalcTrajPosTriHldNegLin,
                        CalcTrajPosLinHldPosTri, s);
  }

  inline static SyncDOFsMotionState D67(const SyncDOFsMotionState& s) {
    // if a->amax->0 and then a->amin->0, is v<=vtrgt?
    return MakeDecision(
        67, AUpToAMaxDownToZero | ADownToAMinUpToZero | CmpVLteVTrgt, &D51,
        &D48, s);
  }

  inline static SyncDOFsMotionState D68(const SyncDOFsMotionState& s) {
    // if a->0, then v->vtrgt(NegTrap) and then a(=0)->hold so that t=tsync, is
    // p>=ptrgt?
    return MakeDecision(68,
                        ADownToZero | VDownToVTrgtNegTrap | HoldToTSync |
                            CmpPGtePTrgtVMinEpsilon,
                        CalcTrajNegLinAdowntoZero | FlipState | &D12, &D69, s);
  }

  inline static SyncDOFsMotionState D69(const SyncDOFsMotionState& s) {
    // if a->0, then a(=0)->hold and then v->vtrgt(NegTrap) so that t=tsync, is
    // p>=(epsilon)ptrgt?
    return MakeDecision(69,
                        ADownToZero | HoldToTSyncAndThen(VDownToVTrgtNegTrap) |
                            CmpPGteEpsilonPTrgt,
                        CalcTrajNegLinAdowntoZero | FlipState | &D28, &D39, s);
  }

  inline static SyncDOFsMotionState D70(const SyncDOFsMotionState& s) {
    // if the NegLinHldPosTri profile was applied, would a valid solution be
    // possible?  If so keep it.
    return TryApplyProfileFunction(70, CalcTrajNegLinHldPosTri,
                                   CalcTrajNegLinHldPosTrap, s);
  }

  inline static SyncDOFsMotionState D71(const SyncDOFsMotionState& s) {
    // if the NegLinHldPosTri profile was applied, would a valid solution be
    // possible?  If so keep it.
    return TryApplyProfileFunction(71, CalcTrajPosLinHldPosTri,
                                   CalcTrajPosLinHldPosTrap, s);
  }

  inline static SyncDOFsMotionState D72(const SyncDOFsMotionState& s) {
    return MakeDecision(72, AUpToAMaxDownToZero | CmpVLteVMin,
                        CalcTrajPosTrapVuptoVMin | &D9,
                        CalcTrajPosTriVuptoVMin | &D9, s);
  }
};

}  // namespace

bool SynchronizeDOFs(const Inputs& inputs, const PositionFlags& flags,
                     PositionOutputs& outputs) {
  SyncDOFsFlags step2_flags = {
      .sync_behavior = flags.synchronization_behavior,
      .use_positional_checks = true,
      .get_into_boundaries_fast =
          flags.behavior_if_initial_state_breaches_constraints ==
          PositionFlags::BehaviorIfInitialStateBreachesConstraints::
              kGetIntoBoundariesFast};

  return SynchronizeDOFs(inputs, step2_flags, &Decisions::D1, outputs);
}

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic
