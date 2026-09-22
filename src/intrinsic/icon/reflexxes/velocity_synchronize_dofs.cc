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

#include <algorithm>
#include <array>

#include "intrinsic/icon/reflexxes/inputs.h"
#include "intrinsic/icon/reflexxes/internal/constants.h"
#include "intrinsic/icon/reflexxes/internal/decision_tree_utility_functions.h"
#include "intrinsic/icon/reflexxes/internal/intermediate_profile_funcs.h"
#include "intrinsic/icon/reflexxes/internal/math.h"
#include "intrinsic/icon/reflexxes/internal/motion_state.h"
#include "intrinsic/icon/reflexxes/internal/motion_state_transforms.h"
#include "intrinsic/icon/reflexxes/internal/motion_state_with_polys.h"
#include "intrinsic/icon/reflexxes/internal/polynomial_solvers.h"
#include "intrinsic/icon/reflexxes/internal/synchronize_dofs_base.h"
#include "intrinsic/icon/reflexxes/profile.h"
#include "intrinsic/icon/reflexxes/velocity_flags.h"
#include "intrinsic/icon/reflexxes/velocity_outputs.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {
namespace {

constexpr auto FinishBlueProfile =
    Curry<SyncDOFsMotionState, Profile, IsMotionState>(
        [](const SyncDOFsMotionState& s_orig, Profile profile, const auto& s) {
          // Check for failure
          if ((fabs(GetV(s) - GetVTrgt(s)) >
               (fabs(kRelVelocityErrorTolerance *
                     (GetV(s_orig) - GetVTrgt(s))) +
                kAbsVelocityErrorTolerance))) {
            if (GetTSync(s) < kMaxAllowedSyncTimeSeconds) {
              return SetFailure(s_orig);
            }
          }

          // Finish off the polynomials by leaving one that goes on forever
          double flip_mult = IsFlipped(s) ? -1.0 : 1.0;
          GetPolynomials(s)->TrimSegmentCount(GetPolynomialIndex(s));
          GetPolynomials(s)->AddSegment({.position = flip_mult * GetP(s),
                                         .velocity = flip_mult * GetVTrgt(s),
                                         .acceleration = 0.,
                                         .jerk = 0.,
                                         .start_time = GetT(s),
                                         .end_time = kInfinity});
          return SetPolynomialIndex(GetPolynomials(s)->GetSegmentCount(), s) |
                 SetProfile(profile);
        });

// NegLinHldNegLin
SyncDOFsMotionState CalcTrajNegLinHldNegLin(const SyncDOFsMotionState& s) {
  return GeneratePolynomials(
             ADownToAHldThenHoldToVTrgtTsyncAndThen(ADownToZero), s) |
         FinishBlueProfile(s, Profile::kNegLinHldNegLin);
}

// PosLinHldNegLin
SyncDOFsMotionState CalcTrajPosLinHldNegLin(const SyncDOFsMotionState& s) {
  // Calculate the hold acceleration.
  std::array<double, 3> coefficients;
  coefficients[2] = 0.5 * GetInvJMaxSubInvJMin(s);
  coefficients[1] = -(GetTSync(s) - GetT(s)) - GetA(s) / GetJMax(s);
  coefficients[0] =
      0.5 * Power2(GetA(s)) / GetJMax(s) + (GetVTrgt(s) - GetV(s));

  double hold_acceleration = GetA(s);
  for (double root : SafeCalculatePolynomialRoots(coefficients)) {
    if ((root - GetA(s)) > -kAbsAccelerationErrorTolerance &&
        (root - GetAMax(s)) < kAbsAccelerationErrorTolerance &&
        (GetTSync(s) - GetT(s) - root * GetInvJMaxSubInvJMin(s) +
         GetA(s) / GetJMax(s)) > -kStep1TimeEpsilon) {
      hold_acceleration = root;
      break;
    }
  }
  hold_acceleration =
      ClipToRange(hold_acceleration, std::max(GetA(s), 0.), GetAMax(s));

  return GeneratePolynomials(
             AUpToAPrime(hold_acceleration) | HoldToTSyncAndThen(ADownToZero),
             s) |
         FinishBlueProfile(s, Profile::kPosLinHldNegLin);
}

// The decision tree.
// This is done as a struct to allow the decisions to be listed in order
// starting from the first decision.  Bundling the decisions into a struct
// allows us to refer to later decisions before they are declared, whereas a
// namespace wouldn't allow that.
struct Decisions {
  static inline SyncDOFsMotionState D1(const SyncDOFsMotionState& s) {
    return MakeDecision(1, CmpAGteZero, &D2, FlipState | &D2, s);
  }

  static inline SyncDOFsMotionState D2(const SyncDOFsMotionState& s) {
    return MakeDecision(2, CmpALteAMax, &D3, CalcTrajNegLinAdowntoAMax | &D3,
                        s);
  }

  static inline SyncDOFsMotionState D3(const SyncDOFsMotionState& s) {
    return MakeDecision(3, ADownToZero | CmpVLteVTrgt, &D4,
                        CalcTrajNegLinAdowntoZero | FlipState | &D4, s);
  }

  static inline SyncDOFsMotionState D4(const SyncDOFsMotionState& s) {
    return MakeDecision(4, HoldToTSyncAndThen(ADownToZero) | CmpVGtVTrgt,
                        CalcTrajNegLinHldNegLin, CalcTrajPosLinHldNegLin, s);
  }
};

}  // namespace

bool SynchronizeDOFs(const Inputs& inputs, const VelocityFlags& flags,
                     VelocityOutputs& outputs) {
  SyncDOFsFlags step2_flags = {.sync_behavior = flags.synchronization_behavior,
                               .use_positional_checks = false,
                               .get_into_boundaries_fast = false};
  return SynchronizeDOFs(inputs, step2_flags, &Decisions::D1, outputs);
}

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic
