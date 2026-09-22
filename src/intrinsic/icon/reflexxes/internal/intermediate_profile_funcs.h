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

#ifndef INTRINSIC_ICON_REFLEXXES_INTERNAL_INTERMEDIATE_PROFILE_FUNCS_H_
#define INTRINSIC_ICON_REFLEXXES_INTERNAL_INTERMEDIATE_PROFILE_FUNCS_H_

#include "intrinsic/icon/reflexxes/internal/motion_state.h"
#include "intrinsic/icon/reflexxes/internal/motion_state_transforms.h"
#include "intrinsic/icon/reflexxes/internal/motion_state_with_polys.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

constexpr auto CalcExecTimeNegLinAdowntoAMax =
    ADownToAMax | AppendProfileTrace(Profile::kNegLinADownToAMax);

constexpr auto CalcTrajNegLinAdowntoAMax =
    GeneratePolynomials(CalcExecTimeNegLinAdowntoAMax);

constexpr auto CalcExecTimeNegLinAdowntoZero =
    ADownToZero | AppendProfileTrace(Profile::kNegLinADownToZero);

constexpr auto CalcTrajNegLinAdowntoZero =
    GeneratePolynomials(CalcExecTimeNegLinAdowntoZero);

constexpr auto CalcExecTimePosLinVuptoVMin =
    VUpToVMinPosLin | AppendProfileTrace(Profile::kPosLinVUpToVMin);

constexpr auto CalcTrajPosLinVuptoVMin =
    GeneratePolynomials(CalcExecTimePosLinVuptoVMin);

constexpr auto CalcExecTimePosLinHldVuptoVMin =
    AUpToAMax | HoldToVMin | AppendProfileTrace(Profile::kPosLinHldVUpToVMin);

constexpr auto CalcTrajPosLinHldVuptoVMin =
    GeneratePolynomials(CalcExecTimePosLinHldVuptoVMin);

// v->vmin (PosLinNegLin) so that v=vmax after a->0
constexpr auto CalcExecTimePosLinNegLinVuptoVMin =
    AUpToAPrime(CalcPeakAccelVMaxPosTri) | VDownToVMinNegLin |
    AppendProfileTrace(Profile::kPosLinNegLinVUpToVMin);

constexpr auto CalcTrajPosLinNegLinVuptoVMin =
    GeneratePolynomials(CalcExecTimePosLinNegLinVuptoVMin);

// PosLinHoldNegLin v increased to vmin
// (PoslinHoldNegLin == PosTrap)
constexpr auto CalcExecTimePosLinHldNegLinVuptoVMin =
    VUpToVMaxPosTrap | VDownToVMinNegLin |
    AppendProfileTrace(Profile::kPosLinHldNegLinVUpToVMin);

// A subtle difference from the Green Profile version in that we do not actually
// decrease a to 0, but account for it later.
constexpr auto CalcTrajPosLinHldNegLinVuptoVMin =
    GeneratePolynomials(AUpToAMax | HoldToVMaxAccountingFor(ADownToZero) |
                        VDownToVMinNegLin) |
    AppendProfileTrace(Profile::kPosLinHldNegLinVUpToVMin);

constexpr auto CalcExecTimePosTrapVuptoVMin =
    VUpToVMinPosTrap | AppendProfileTrace(Profile::kPosTrapVUpToVMin);

constexpr auto CalcTrajPosTrapVuptoVMin =
    GeneratePolynomials(CalcExecTimePosTrapVuptoVMin);

constexpr auto CalcExecTimePosTriVuptoVMin =
    VUpToVMinPosTri | AppendProfileTrace(Profile::kPosTriVUpToVMin);

constexpr auto CalcTrajPosTriVuptoVMin =
    GeneratePolynomials(CalcExecTimePosTriVuptoVMin);

constexpr auto CalcExecTimePosTrapVuptoVMax =
    VUpToVMaxPosTrap | AppendProfileTrace(Profile::kPosTrapVUpToVMax);

constexpr auto CalcExecTimePosTriVuptoVMax =
    VUpToVMaxPosTri | AppendProfileTrace(Profile::kPosTriVUpToVMax);

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_INTERNAL_INTERMEDIATE_PROFILE_FUNCS_H_
