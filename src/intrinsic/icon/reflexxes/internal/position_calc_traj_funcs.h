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

#ifndef INTRINSIC_ICON_REFLEXXES_INTERNAL_POSITION_CALC_TRAJ_FUNCS_H_
#define INTRINSIC_ICON_REFLEXXES_INTERNAL_POSITION_CALC_TRAJ_FUNCS_H_

#include <cmath>
#include <cstring>
#include <vector>

#include "intrinsic/icon/reflexxes/inputs.h"
#include "intrinsic/icon/reflexxes/internal/constants.h"
#include "intrinsic/icon/reflexxes/internal/motion_state_transforms.h"
#include "intrinsic/icon/reflexxes/internal/motion_state_with_polys.h"
#include "intrinsic/icon/reflexxes/internal/synchronize_dofs_base.h"
#include "intrinsic/icon/reflexxes/motion_polynomials.h"
#include "intrinsic/icon/reflexxes/position_inputs.h"
#include "intrinsic/icon/reflexxes/position_outputs.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

SyncDOFsMotionState CalcTrajPosTrapZeroNegTrap(const SyncDOFsMotionState& s);

SyncDOFsMotionState CalcTrajPosTrapZeroNegTri(const SyncDOFsMotionState& s);

SyncDOFsMotionState CalcTrajPosTriZeroNegTrap(const SyncDOFsMotionState& s);

SyncDOFsMotionState CalcTrajPosTriZeroNegTri(const SyncDOFsMotionState& s);

SyncDOFsMotionState CalcTrajPosTrapZeroPosTrap(const SyncDOFsMotionState& s);

SyncDOFsMotionState CalcTrajPosTrapZeroPosTri(const SyncDOFsMotionState& s);

SyncDOFsMotionState CalcTrajPosTriZeroPosTrap(const SyncDOFsMotionState& s);

SyncDOFsMotionState CalcTrajPosTriZeroPosTri(const SyncDOFsMotionState& s);

SyncDOFsMotionState CalcTrajPosLinHldPosTrap(const SyncDOFsMotionState& s);

SyncDOFsMotionState CalcTrajPosLinHldPosTri(const SyncDOFsMotionState& s);

SyncDOFsMotionState CalcTrajNegLinHldPosTrap(const SyncDOFsMotionState& s);

SyncDOFsMotionState CalcTrajNegLinHldPosTri(const SyncDOFsMotionState& s);

SyncDOFsMotionState CalcTrajPosTrapHldNegLin(const SyncDOFsMotionState& s);

SyncDOFsMotionState CalcTrajPosTriHldNegLin(const SyncDOFsMotionState& s);

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_INTERNAL_POSITION_CALC_TRAJ_FUNCS_H_
