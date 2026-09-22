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

#ifndef INTRINSIC_ICON_REFLEXXES_PROFILE_H_
#define INTRINSIC_ICON_REFLEXXES_PROFILE_H_

#include <string>

#include "absl/strings/str_cat.h"
namespace intrinsic {
namespace reflexxes {

// Describes the choice of a motion profile by the algorithm.
enum class Profile : int8_t {
  // Not set or error.
  kUndefined,

  // The set of "final" (red/blue) profiles.
  kPosTrapNegTrap,
  kPosTrapNegTri,
  kPosTriNegTrap,
  kPosTriNegTri,
  kPosTrapZeroNegTrap,
  kPosTrapZeroNegTri,
  kPosTriZeroNegTrap,
  kPosTriZeroNegTri,
  kNegLinPosTrap,
  kNegLinPosTri,
  kPosTrapZeroPosTrap,
  kPosTrapZeroPosTri,
  kPosTriZeroPosTrap,
  kPosTriZeroPosTri,
  kPosLinHldPosTrap,
  kPosLinHldPosTri,
  kNegLinHldPosTrap,
  kNegLinHldPosTri,
  kPosTrapHldNegLin,
  kPosTriHldNegLin,
  kPosTrap,
  kPosTri,
  kNegLinHldNegLin,
  kPosLinHldNegLin,

  // The following are only used as intermediate "green"/"pink" profiles.
  kNegLinADownToAMax,
  kNegLinADownToZero,
  kPosLinVUpToVMin,
  kPosLinHldVUpToVMin,
  kPosLinNegLinVUpToVMin,
  kPosLinHldNegLinVUpToVMin,
  kPosTrapVUpToVMin,
  kPosTrapVUpToVMax,
  kPosTriVUpToVMin,
  kPosTriVUpToVMax,

  // Also only used as intermediate, indicates the profile was flipped.
  kFlip,

  // Indicates the profile was overridden during phase sync.
  kOverridden,
};

constexpr const char* GetProfileString(const Profile profile) {
  switch (profile) {
    case Profile::kUndefined:
      return "Undefined";
    case Profile::kPosTrapNegTrap:
      return "PosTrapNegTrap";
    case Profile::kPosTrapNegTri:
      return "PosTrapNegTri";
    case Profile::kPosTriNegTrap:
      return "PosTriNegTrap";
    case Profile::kPosTriNegTri:
      return "PosTriNegTri";
    case Profile::kPosTrapZeroNegTrap:
      return "PosTrapZeroNegTrap";
    case Profile::kPosTrapZeroNegTri:
      return "PosTrapZeroNegTri";
    case Profile::kPosTriZeroNegTrap:
      return "PosTriZeroNegTrap";
    case Profile::kPosTriZeroNegTri:
      return "PosTriZeroNegTri";
    case Profile::kNegLinPosTrap:
      return "NegLinPosTrap";
    case Profile::kNegLinPosTri:
      return "NegLinPosTri";
    case Profile::kPosTrapZeroPosTrap:
      return "PosTrapZeroPosTrap";
    case Profile::kPosTrapZeroPosTri:
      return "PosTrapZeroPosTri";
    case Profile::kPosTriZeroPosTrap:
      return "PosTriZeroPosTrap";
    case Profile::kPosTriZeroPosTri:
      return "PosTriZeroPosTri";
    case Profile::kPosLinHldPosTrap:
      return "PosLinHldPosTrap";
    case Profile::kPosLinHldPosTri:
      return "PosLinHldPosTri";
    case Profile::kNegLinHldPosTrap:
      return "NegLinHldPosTrap";
    case Profile::kNegLinHldPosTri:
      return "NegLinHldPosTri";
    case Profile::kPosTrapHldNegLin:
      return "PosTrapHldNegLin";
    case Profile::kPosTriHldNegLin:
      return "PosTriHldNegLin";
    case Profile::kPosTrap:
      return "PosTrap";
    case Profile::kPosTri:
      return "PosTri";
    case Profile::kNegLinHldNegLin:
      return "NegLinHldNegLin";
    case Profile::kPosLinHldNegLin:
      return "PosLinHldNegLin";
    case Profile::kNegLinADownToAMax:
      return "NegLinADownToAMax";
    case Profile::kNegLinADownToZero:
      return "NegLinADownToZero";
    case Profile::kPosLinVUpToVMin:
      return "PosLinVUpToVMin";
    case Profile::kPosLinHldVUpToVMin:
      return "PosLinHldVupToVMin";
    case Profile::kPosLinNegLinVUpToVMin:
      return "PosLinNegLinVUpToVMin";
    case Profile::kPosLinHldNegLinVUpToVMin:
      return "PosLinHldNegLinVUpToVMin";
    case Profile::kPosTrapVUpToVMin:
      return "PosTrapVUpToVMin";
    case Profile::kPosTrapVUpToVMax:
      return "PosTrapVUpToVMax";
    case Profile::kPosTriVUpToVMin:
      return "PosTriVUpToVMin";
    case Profile::kPosTriVUpToVMax:
      return "PosTriVUpToVMax";
    case Profile::kFlip:
      return "Flip";
    case Profile::kOverridden:
      return "Overridden";
  }
}

}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_PROFILE_H_
