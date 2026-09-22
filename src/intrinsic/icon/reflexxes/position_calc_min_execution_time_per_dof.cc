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

// This file defines the decision tree for the CalcMinExecutionTimePerDOF
// function, otherwise known as "Position Step1A" in the Reflexxes literature.
// See
// https://github.com/intrinsic-ai/intrinsic-core/blob/main/intrinsic_control/intrinsic/icon/reflexxes/g3doc/step1a_decision_tree.pdf
// for a detailed description of what's happening in this file.

#include "intrinsic/icon/reflexxes/inputs.h"
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
#include "intrinsic/util/functional.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

namespace {

Range UpdateLimitsForNegLinPosTri(const bool get_min_root,
                                  const Range& param_min_max,
                                  const MotionState& s) {
  auto p_func = GetPositionErrorForNegLinPosTriFunc(s);
  auto dp_func = GetDerivativePositionErrorForNegLinPosTri(s);
  auto ddp_func = GetSecondDerivativePositionErrorForNegLinPosTri(s);
  Range result = param_min_max;
  double position_error_at_param_min = p_func(result.min, nullptr, nullptr);
  double position_error_at_param_max = p_func(result.max, nullptr, nullptr);

  double derivative_position_error_at_param_min =
      dp_func(result.min, nullptr, nullptr);
  double derivative_position_error_at_param_max =
      dp_func(result.max, nullptr, nullptr);

  // check whether derivative changes sign
  if (GetSign(derivative_position_error_at_param_min) !=
      GetSign(derivative_position_error_at_param_max)) {
    // use bisection method to find the extremum
    double param_extremum = GetRootABKMethod(dp_func, result);

    // position error at extremum
    double position_error_at_param_extremum =
        p_func(param_extremum, nullptr, nullptr);

    if (get_min_root) {
      if (GetSign(position_error_at_param_extremum) !=
              GetSign(position_error_at_param_min) &&
          param_extremum > result.min) {
        result.max = param_extremum;
        return result;
      }
    } else if (GetSign(position_error_at_param_extremum) !=
                   GetSign(position_error_at_param_max) &&
               param_extremum < result.max) {
      result.min = param_extremum;
      return result;
    }
  } else {
    double second_derivative_position_error_at_param_min =
        ddp_func(result.min, nullptr, nullptr);
    double second_derivative_position_error_at_param_max =
        ddp_func(result.max, nullptr, nullptr);

    // check whether the second derivative changes sign
    if (GetSign(second_derivative_position_error_at_param_min) !=
        GetSign(second_derivative_position_error_at_param_max)) {
      // use bisection method to find the extremum of derivative position
      // error
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
          if ((GetSign(position_error_at_param_lower_extremum) !=
               GetSign(position_error_at_param_upper_extremum)) &&
              (GetSign(position_error_at_param_lower_extremum) !=
               GetSign(position_error_at_param_min)) &&
              param_lower_extremum > result.min) {
            result.max = param_lower_extremum;
            return result;
          }
        } else {
          if ((GetSign(position_error_at_param_lower_extremum) !=
               GetSign(position_error_at_param_upper_extremum)) &&
              (GetSign(position_error_at_param_max) !=
               GetSign(position_error_at_param_upper_extremum)) &&
              param_upper_extremum < result.max) {
            result.min = param_upper_extremum;
            return result;
          }
        }
      }
    }
  }
  return result;
}

MotionState CalcExecTimeNegLinPosTri(const MotionState& s) {
  auto param_min_max = UpdateLimitsForNegLinPosTri(
      false, GetNegLinEndAccelerationLimitsForNegLinPosTri(s), s);

  // return execution time
  return RedProfile(
      Profile::kNegLinPosTri, GetPositionErrorForNegLinPosTriFunc(s),
      GetDerivativePositionErrorForNegLinPosTri(s), false, param_min_max, s);
}

inline MotionState CalcExecTimePosTriNegTri(const MotionState& s) {
  return RedProfile(Profile::kPosTriNegTri, GetPositionErrorForPosTriNegTri(s),
                    GetDerivativePositionErrorForPosTriNegTri(s), true,
                    GetPosTriPeakAccelerationLimitsForPosTriNegTri(s), s);
}

constexpr auto ApplyCalcExecTimeIfElse = functional::CurryN<3>(
    [](auto red_profile_func1, auto red_profile_func2, const MotionState& s) {
      MotionState s1 = red_profile_func1(s);
      if (IsSuccess(s1)) {
        return s1;
      }
      return red_profile_func2(s);
    });

// Flips the state, and all sets our special flag to trace whether the input is
// flipped.
constexpr auto FlipStateAndTrace = [](const MotionState& s) {
  return FlipState(s) |
         SetStep1AAppliedProfileFlipped(!GetStep1AAppliedProfileFlipped(s));
};

// The decision tree.
// This is done as a struct to allow the decisions to be listed in order
// starting from the first decision.  Bundling the decisions into a struct
// allows us to refer to later decisions before they are declared, whereas a
// namespace wouldn't allow that.
struct Decisions {
  inline static MotionState D1(const MotionState& s) {
    // Note the explicit non-use of FlipStateAndTrace.
    return MakeDecision(1, CmpAGteZero, &D2, FlipState | &D2, s);
  }

  inline static MotionState D2(const MotionState& s) {
    return MakeDecision(2, CmpALteAMax, &D3,
                        CalcExecTimeNegLinAdowntoAMax | &D3, s);
  }

  inline static MotionState D3(const MotionState& s) {
    // Note the explicit non-use of FlipStateAndTrace.
    return MakeDecision(3, ADownToZero | CmpVLteVMax, &D4,
                        CalcExecTimeNegLinAdowntoZero | FlipState |
                            IfElse(GetIntoBoundariesFast(s), &D5, &D37),
                        s);
  }

  inline static MotionState D4(const MotionState& s) {
    return MakeDecision(4, CmpVGteVMin, &D9,
                        IfElse(GetIntoBoundariesFast(s), &D5, &D37), s);
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
    // if v->vmin (PosLinHld) and then a->0, is v>vmax?
    // Just start after the PosLinHld, which is v=vmin, a=amax
    return MakeDecision(
        7, SetV(GetVMin(s), s) | SetA(GetAMax(s)) | ADownToZero | CmpVGtVMax,
        &D8, CalcExecTimePosLinHldVuptoVMin | &D9, s);
  }

  inline static MotionState D8(const MotionState& s) {
    return MakeDecision(8, AUpToAMaxDownToZero | CmpVGtVMax,
                        CalcExecTimePosLinNegLinVuptoVMin | &D9,
                        CalcExecTimePosLinHldNegLinVuptoVMin | &D9, s);
  }

  inline static MotionState D9(const MotionState& s) {
    return MakeDecision(9, ADownToZero | CmpVLteVTrgt, &D10, &D33, s);
  }

  inline static MotionState D10(const MotionState& s) {
    return MakeDecision(10, AUpToAMaxDownToZero | CmpVLteVTrgt, &D11, &D31, s);
  }

  inline static MotionState D11(const MotionState& s) {
    return MakeDecision(11, VUpToVTrgtPosTrap | CmpPLtePTrgt, &D12, &D16, s);
  }

  inline static MotionState D12(const MotionState& s) {
    return MakeDecision(12,
                        VUpToVMaxPosTrap | ADownToAMinUpToZero | CmpVLteVTrgt,
                        &D13, &D14, s);
  }

  inline static MotionState D13(const MotionState& s) {
    return MakeDecision(13,
                        VUpToVMaxPosTrap | VDownToVTrgtNegTri | CmpPLtePTrgt,
                        CalcExecTimePosTrapZeroNegTri,
                        CalcExecTimePosTrapNegTri(false, true, 0), s);
  }

  inline static MotionState D14(const MotionState& s) {
    return MakeDecision(
        14,
        AUpToAMax | HoldToVTrgtAndThen(ADownToZero | ADownToAMinUpToZero) |
            CmpPLtePTrgt,
        &D15, CalcExecTimePosTrapNegTri(false, true, 0), s);
  }

  inline static MotionState D15(const MotionState& s) {
    return MakeDecision(15,
                        VUpToVMaxPosTrap | VDownToVTrgtNegTrap | CmpPLtePTrgt,
                        CalcExecTimePosTrapZeroNegTrap,
                        CalcExecTimePosTrapNegTrap(true, true, 0), s);
  }

  inline static MotionState D16(const MotionState& s) {
    return MakeDecision(16, ADownToZero | AUpToAMaxDownToZero | CmpVLteVTrgt,
                        &D17, &D20, s);
  }

  inline static MotionState D17(const MotionState& s) {
    return MakeDecision(
        17, ADownToZero | VUpToVTrgtPosTrap | CmpPLtePTrgt,
        CalcExecTimeNegLinPosTrap(false, true, 0),
        CalcExecTimeNegLinAdowntoZero | FlipStateAndTrace | &D18, s);
  }

  inline static MotionState D18(const MotionState& s) {
    return MakeDecision(18, AUpToAMaxDownToZero | CmpVLteVMax, &D19, &D30, s);
  }

  inline static MotionState D19(const MotionState& s) {
    return MakeDecision(
        19, AUpToAMaxDownToZero | VDownToVTrgtNegTrap | CmpPLtePTrgt, &D15,
        CalcExecTimePosTriNegTrap(true, true, 0), s);
  }

  inline static MotionState D20(const MotionState& s) {
    return MakeDecision(
        20, ADownToZero | VUpToVTrgtPosTri | CmpPLtePTrgt, &D21,
        CalcExecTimeNegLinAdowntoZero | FlipStateAndTrace | &D22, s);
  }

  inline static MotionState D21(const MotionState& s) {
    // if the NegLinPosTrap profile was applied, would a valid solution for tmin
    // possible?
    return MakeDecision(
        21, CalcExecTimeNegLinPosTrap(false, false, 0) | IsSuccess,
        ApplyCalcExecTimeIfElse(CalcExecTimeNegLinPosTrap(false, true, 0),
                                &CalcExecTimeNegLinPosTri),
        ApplyCalcExecTimeIfElse(&CalcExecTimeNegLinPosTri,
                                CalcExecTimeNegLinPosTrap(false, true, 0)),
        s);
  }

  inline static MotionState D22(const MotionState& s) {
    return MakeDecision(
        22, AUpToAMaxDownToZero | ADownToAMinUpToZero | CmpVLteVTrgt, &D23,
        &D26, s);
  }

  inline static MotionState D23(const MotionState& s) {
    return MakeDecision(23, AUpToAMaxDownToZero | CmpVLteVMax, &D36, &D28, s);
  }

  inline static MotionState D24(const MotionState& s) {
    return MakeDecision(24,
                        AUpToAMaxDownToZero | VDownToVTrgtNegTri | CmpPLtePTrgt,
                        &D12, &CalcExecTimePosTriNegTri, s);
  }

  inline static MotionState D25(const MotionState& s) {
    return MakeDecision(25, VUpToVMaxPosTri | VDownToVTrgtNegTri | CmpPLtePTrgt,
                        CalcExecTimePosTriZeroNegTri, &CalcExecTimePosTriNegTri,
                        s);
  }

  inline static MotionState D26(const MotionState& s) {
    return MakeDecision(26, AUpToAMaxDownToZero | CmpVLteVMax, &D27, &D28, s);
  }

  inline static MotionState D27(const MotionState& s) {
    // if a->apeak->0 and then a->amin->0, so that v=vtrgt, is p<=ptrgt?
    return MakeDecision(
        27, VToVTrgtPosTriAndThen(ADownToAMinUpToZero, s) | CmpPLtePTrgt, &D19,
        &CalcExecTimePosTriNegTri, s);
  }

  inline static MotionState D28(const MotionState& s) {
    return MakeDecision(28,
                        VUpToVMaxPosTri | ADownToAMinUpToZero | CmpVLteVTrgt,
                        &D25, &D29, s);
  }

  inline static MotionState D29(const MotionState& s) {
    // if a->apeak->0 and then a->amin->0, so that v=vtrgt, is p<=ptrgt?
    return MakeDecision(
        29, VToVTrgtPosTriAndThen(ADownToAMinUpToZero, s) | CmpPLtePTrgt, &D30,
        &CalcExecTimePosTriNegTri, s);
  }

  inline static MotionState D30(const MotionState& s) {
    return MakeDecision(30,
                        VUpToVMaxPosTri | VDownToVTrgtNegTrap | CmpPLtePTrgt,
                        CalcExecTimePosTriZeroNegTrap,
                        CalcExecTimePosTriNegTrap(true, true, 0), s);
  }

  inline static MotionState D31(const MotionState& s) {
    return MakeDecision(31, VUpToVTrgtPosTri | CmpPLtePTrgt, &D23, &D32, s);
  }

  inline static MotionState D32(const MotionState& s) {
    return MakeDecision(
        32, ADownToZero | VUpToVTrgtPosTri | CmpPLtePTrgt,
        &CalcExecTimeNegLinPosTri,
        CalcExecTimeNegLinAdowntoZero | FlipStateAndTrace | &D22, s);
  }

  inline static MotionState D33(const MotionState& s) {
    return MakeDecision(33, ADownToZero | ADownToAMinUpToZero | CmpVLteVTrgt,
                        &D34, &D35, s);
  }

  inline static MotionState D34(const MotionState& s) {
    return MakeDecision(
        34, ADownToZero | VDownToVTrgtNegTri | CmpPLtePTrgt, &D22,
        CalcExecTimeNegLinAdowntoZero | FlipStateAndTrace | &D23, s);
  }

  inline static MotionState D35(const MotionState& s) {
    return MakeDecision(
        35, ADownToZero | VDownToVTrgtNegTrap | CmpPLtePTrgt, &D18,
        CalcExecTimeNegLinAdowntoZero | FlipStateAndTrace | &D12, s);
  }

  inline static MotionState D36(const MotionState& s) {
    return MakeDecision(
        36, AUpToAMaxDownToZero | ADownToAMinUpToZero | CmpVLteVTrgt, &D24,
        &D27, s);
  }

  inline static MotionState D37(const MotionState& s) {
    return MakeDecision(37, AUpToAMaxDownToZero | CmpVLteVMin,
                        CalcExecTimePosTrapVuptoVMin | &D9,
                        CalcExecTimePosTriVuptoVMin | &D9, s);
  }
};

}  // namespace

bool CalcMinExecutionTimePerDOF(const PositionInputs& inputs,
                                const PositionFlags& flags,
                                PositionOutputs& outputs,
                                MaxDOFFixedVector<bool>& step1a_flipped) {
  for (auto [dof_input, dof_output] :
       Zip(inputs.GetDOFs(), outputs.GetDOFs())) {
    // Check if the dof is available for use and if the constraints are valid
    if (!dof_input.selected) {
      dof_output.step1a.result = true;
      continue;
    }

    MotionState s = RunPositionStep1DecisionTree(
        Decisions::D1, dof_input, flags, dof_output, dof_output.step1a);

    if (IsFailure(s)) {
      // an error occurred so quit deciding for the rest of the dofs
      // and notify failure
      return false;
    }

    dof_output.min_execution_time = GetT(s);
    WriteStateToSubStepOutput(s, dof_output.step1a);

    // Save off whether or not we flipped.
    step1a_flipped.at(dof_input.index) = GetStep1AAppliedProfileFlipped(s);
  }

  return true;
}

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic
