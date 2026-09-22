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

// This file defines the decision tree for the CalcInoperativeZoneEnd
// function, otherwise known as "Position Step1C" in the Reflexxes literature.
// See
// https://github.com/intrinsic-ai/intrinsic-core/blob/main/intrinsic_control/intrinsic/icon/reflexxes/g3doc/step1c_decision_tree.pdf
// for a detailed description of what's happening in this file.

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

// Template for updating the parameters limits for Step 1C
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
          if (param_extremum < param_min_max.max) {
            return {param_extremum, param_min_max.max};
          }
        } else {
          if (param_extremum > param_min_max.min) {
            return {param_min_max.min, param_extremum};
          }
        }
      }
    }
  }
  return param_min_max;
}

constexpr auto UpdateLimitsForPosTri = Curry(UpdateLimitsForPosTriFunc);

MotionState CalcExecTimePosTriNegTriFunc(
    const double step1b_red_profile_execution_time, const MotionState& s) {
  return CalcExecTimePosTriNegTri(UpdateLimitsForPosTri,
                                  step1b_red_profile_execution_time, s);
}
constexpr auto CalcExecTimePosTriNegTri = Curry(CalcExecTimePosTriNegTriFunc);

double GetCalcExecTimeExecutionTime(const MotionState& s) {
  return GetDOFOutput(s).inoperative_begin_execution_time - GetT(s);
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
                            IfElse(GetIntoBoundariesFast(s), &D5, &D24),
                        s);
  }

  inline static MotionState D4(const MotionState& s) {
    return MakeDecision(4, CmpVGteVMin, &D9,
                        IfElse(GetIntoBoundariesFast(s), &D5, &D24), s);
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
    // Note that we set the conditions _after_ PosLinHld, assuming a->amax
    return MakeDecision(
        7, SetV(GetVMin(s), s) | SetA(GetAMax(s)) | ADownToZero | CmpVGtVMax,
        &D8, CalcExecTimePosLinHldVuptoVMin | &D9, s);
  }

  inline static MotionState D8(const MotionState& s) {
    return MakeDecision(8, AUpToAMaxDownToZero(s) | CmpVGtVMax,
                        CalcExecTimePosLinNegLinVuptoVMin | &D9,
                        CalcExecTimePosLinHldNegLinVuptoVMin | &D9, s);
  }

  inline static MotionState D9(const MotionState& s) {
    return MakeDecision(9, CmpVTrgtLteZero, &D10, FlipState | &D10, s);
  }

  inline static MotionState D10(const MotionState& s) {
    return MakeDecision(
        10, CmpAGteZero, &D11,
        FlipState | CalcExecTimeNegLinAdowntoZero | FlipState | &D12, s);
  }

  inline static MotionState D11(const MotionState& s) {
    // Was the result of CalcMinExecTime either NegLinPosTri or NegLinPosTrap
    // AND Was the result of CalcInoperativeZoneStart either NegLinPosTri or
    // NegLinPosTrap ?
    bool decision =
        ((GetDOFOutput(s).step1a.applied_profile == Profile::kNegLinPosTri ||
          GetDOFOutput(s).step1a.applied_profile == Profile::kNegLinPosTrap) &&
         (GetDOFOutput(s).step1b.applied_profile == Profile::kNegLinPosTri ||
          GetDOFOutput(s).step1b.applied_profile == Profile::kNegLinPosTrap));
    return MakeDecision(11, decision, CalcExecTimeNegLinAdowntoZero | &D12,
                        &D23, s);
  }

  inline static MotionState D12(const MotionState& s) {
    return MakeDecision(12, AUpToAMaxDownToZero | CmpVLteVMax, &D13, &D20, s);
  }

  inline static MotionState D13(const MotionState& s) {
    return MakeDecision(13,
                        VUpToVMaxPosTrap | ADownToAMinUpToZero | CmpVLteVTrgt,
                        &D14, &D16, s);
  }

  inline static MotionState D14(const MotionState& s) {
    return MakeDecision(14,
                        VUpToVMaxPosTrap | VDownToVTrgtNegTri | CmpPLtePTrgt,
                        CalcExecTimePosTrapZeroNegTri, &D15, s);
  }

  inline static MotionState D15(const MotionState& s) {
    // if the PosTrapNegTri profile was applied, would a valid solution
    // for tmin possible?
    return TryApplyProfileFunction(
        15,
        CalcExecTimePosTrapNegTri(true, true, GetCalcExecTimeExecutionTime(s)),
        CalcExecTimePosTriNegTri(GetCalcExecTimeExecutionTime(s)), s);
  }

  inline static MotionState D16(const MotionState& s) {
    return MakeDecision(16,
                        VUpToVMaxPosTrap | VDownToVTrgtNegTrap | CmpPLtePTrgt,
                        CalcExecTimePosTrapZeroNegTrap, &D17, s);
  }

  inline static MotionState D17(const MotionState& s) {
    // if the PosTrapNegTrap profile was applied, would a valid solution
    // for tmin possible? conditions: amin<=a<=0 and 0<=v<=vtrgt
    return TryApplyProfileFunction(
        17,
        CalcExecTimePosTrapNegTrap(false, true,
                                   GetCalcExecTimeExecutionTime(s)),
        &D18, s);
  }

  inline static MotionState D18(const MotionState& s) {
    return MakeDecision(18,
                        VUpToVMaxPosTrap | ADownToAMinUpToZero | CmpVLteVTrgt,
                        &D15, &D19, s);
  }

  inline static MotionState D19(const MotionState& s) {
    // if the PosTriNegTrap profile was applied, would a valid solution
    // for tmin possible?
    return TryApplyProfileFunction(
        19,
        CalcExecTimePosTriNegTrap(false, true, GetCalcExecTimeExecutionTime(s)),
        CalcExecTimePosTriNegTri(GetCalcExecTimeExecutionTime(s)), s);
  }

  inline static MotionState D20(const MotionState& s) {
    return MakeDecision(20,
                        VUpToVMaxPosTri | ADownToAMinUpToZero | CmpVLteVTrgt,
                        &D21, &D22, s);
  }

  inline static MotionState D21(const MotionState& s) {
    return MakeDecision(
        21, VUpToVMaxPosTri | VDownToVTrgtNegTri | CmpPLtePTrgt,
        CalcExecTimePosTriZeroNegTri,
        CalcExecTimePosTriNegTri(GetCalcExecTimeExecutionTime(s)), s);
  }

  inline static MotionState D22(const MotionState& s) {
    return MakeDecision(22,
                        VUpToVMaxPosTri | VDownToVTrgtNegTrap | CmpPLtePTrgt,
                        CalcExecTimePosTriZeroNegTrap, &D19, s);
  }

  inline static MotionState D23(const MotionState& s) {
    bool decision =
        (GetDOFOutput(s).step1b.applied_profile == Profile::kNegLinPosTri ||
         GetDOFOutput(s).step1b.applied_profile == Profile::kNegLinPosTrap);
    return MakeDecision(23, decision, CalcExecTimeNegLinAdowntoZero | &D12,
                        &D12, s);
  }

  inline static MotionState D24(const MotionState& s) {
    return MakeDecision(24, AUpToAMaxDownToZero | CmpVLteVMin,
                        CalcExecTimePosTrapVuptoVMin | &D9,
                        CalcExecTimePosTriVuptoVMin | &D9, s);
  }
};

}  // namespace

bool CalcInoperativeZoneEnd(const PositionInputs& inputs,
                            const PositionFlags& flags,
                            PositionOutputs& outputs) {
  for (auto [dof_input, dof_output] :
       Zip(inputs.GetDOFs(), outputs.GetDOFs())) {
    // Check if the dof is available for use and if the constraints are valid.
    // Also check whether an inoperative time interval exists
    // And the solution has not already been computed in
    // CalcInoperativeZoneStart
    if (!dof_input.selected ||
        (dof_output.inoperative_begin_execution_time == kInfinity) ||
        (dof_output.inoperative_end_execution_time != 0.0)) {
      if (dof_output.inoperative_begin_execution_time == kInfinity) {
        dof_output.inoperative_end_execution_time = kInfinity;
      }

      // update the profile error element
      dof_output.step1c.result = true;
      continue;
    }

    MotionState s = RunPositionStep1DecisionTree(
        Decisions::D1, dof_input, flags, dof_output, dof_output.step1c);
    if (IsFailure(s)) {
      // an error occurred so quit deciding for the rest of the dofs
      // and notify failure
      return false;
    }

    dof_output.inoperative_end_execution_time = GetT(s);
    WriteStateToSubStepOutput(s, dof_output.step1c);
  }

  return true;
}

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic
