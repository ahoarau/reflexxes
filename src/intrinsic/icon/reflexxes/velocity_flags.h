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

#ifndef INTRINSIC_ICON_REFLEXXES_VELOCITY_FLAGS_H_
#define INTRINSIC_ICON_REFLEXXES_VELOCITY_FLAGS_H_

#include "intrinsic/icon/reflexxes/flags.h"

namespace intrinsic {
namespace reflexxes {

// Class containing flags to parameterize the execution of the velocity-based
// Online Trajectory Generation algorithm
struct VelocityFlags : public Flags {
  VelocityFlags()
      : Flags(SyncBehavior::kNoSynchronization,
              PositionalLimitsBehavior::kIgnore,
              InvalidConstraintsBehavior::kDeselectDofWithoutErrorMsg,
              InvalidScaleOfInputValuesBehavior::kIgnore) {}

  // Construct from base.
  VelocityFlags(const Flags& flags) : Flags(flags) {}
};

}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_VELOCITY_FLAGS_H_
