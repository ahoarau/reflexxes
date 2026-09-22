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

#ifndef INTRINSIC_ICON_REFLEXXES_POSITION_FLAGS_H_
#define INTRINSIC_ICON_REFLEXXES_POSITION_FLAGS_H_

#include "intrinsic/icon/reflexxes/flags.h"

namespace intrinsic {
namespace reflexxes {

// Class containing flags to parameterize the execution of the position-based
// Online Trajectory Generation algorithm.
struct PositionFlags : public Flags {
  // Enumeration whose values specify the behavior after the final state of
  // motion is reached.
  enum class FinalMotionBehavior {
    // The desired velocity of the target state of motion will be kept at zero
    // acceleration (default).
    kKeepTargetVelocity = 0,
    // After the final state of motion is reached, a new trajectory will be
    // computed, such that the desired state of motion will be reached again
    // (and again, and again, etc.).
    kRecomputeTrajectory = 1,
  };

  // Enumeration whose values specify the behavior if the initial state of
  // motion is out of bounds.
  enum class BehaviorIfInitialStateBreachesConstraints {
    // The system will be brought back into the boundaries with an intermediate
    // acceleration of zero allowing to continue the motion exactly at the
    // maximum velocity. (default).
    kGetIntoBoundariesAtZeroAcceleration = 0,
    // The system will be brought back into the boundaries as fast as possible
    // but with a non-zero intermediate acceleration.
    kGetIntoBoundariesFast = 1,
  };

  PositionFlags()
      : Flags(SyncBehavior::kPhaseSynchronizationIfPossible,
              PositionalLimitsBehavior::kIgnore,
              InvalidConstraintsBehavior::kDeselectDofWithoutErrorMsg,
              InvalidScaleOfInputValuesBehavior::kIgnore) {}

  bool operator==(const PositionFlags& flags) const {
    return ((Flags::operator==(flags)) &&
            (final_motion_behavior == flags.final_motion_behavior) &&
            (behavior_if_initial_state_breaches_constraints ==
             flags.behavior_if_initial_state_breaches_constraints) &&
            (keep_current_velocity_in_case_of_fallback_strategy ==
             flags.keep_current_velocity_in_case_of_fallback_strategy));
  }

  bool operator!=(const PositionFlags& flags) const {
    return !(*this == flags);
  }

  FinalMotionBehavior final_motion_behavior =
      FinalMotionBehavior::kKeepTargetVelocity;

  BehaviorIfInitialStateBreachesConstraints
      behavior_if_initial_state_breaches_constraints =
          BehaviorIfInitialStateBreachesConstraints::
              kGetIntoBoundariesAtZeroAcceleration;

  // If true, PositionInputs::AlternativeTargetVelocityVector will be used in
  // TypeIVRMLPosition::FallBackStrategy().
  bool keep_current_velocity_in_case_of_fallback_strategy = false;
};

}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_POSITION_FLAGS_H_
