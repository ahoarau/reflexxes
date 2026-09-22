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

// A common implementation for the SynchronizeDOFs (aka Step2) modules for both
// Position and Velocity.  These steps are almost identical except for their
// actual decision trees.  All common functionality is grouped into the
// SynchronizeDOFs function, which executes the passed in decision tree along
// with the other tasks like phase synchronization scaling.

#ifndef INTRINSIC_ICON_REFLEXXES_INTERNAL_SYNCHRONIZE_DOFS_BASE_H_
#define INTRINSIC_ICON_REFLEXXES_INTERNAL_SYNCHRONIZE_DOFS_BASE_H_

#include <functional>
#include <vector>

#include "absl/functional/function_ref.h"
#include "intrinsic/icon/reflexxes/inputs.h"
#include "intrinsic/icon/reflexxes/internal/constants.h"
#include "intrinsic/icon/reflexxes/internal/motion_state.h"
#include "intrinsic/icon/reflexxes/motion_polynomials.h"
#include "intrinsic/icon/reflexxes/outputs.h"
#include "intrinsic/icon/reflexxes/position_flags.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

// Subset of Position/Velocity Flags that captures only what Step2 actually
// cares about.
struct SyncDOFsFlags {
  Flags::SyncBehavior sync_behavior =
      Flags::SyncBehavior::kPhaseSynchronizationIfPossible;
  bool use_positional_checks = false;
  bool get_into_boundaries_fast = false;
};

using SyncDOFsMotionState = MotionStateT<>;

// Function that executes a Step2 decision tree
bool SynchronizeDOFs(
    const Inputs& inputs, const SyncDOFsFlags& flags,
    absl::FunctionRef<SyncDOFsMotionState(const SyncDOFsMotionState&)>
        exec_decision_tree,
    Outputs& outputs);

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_INTERNAL_SYNCHRONIZE_DOFS_BASE_H_
