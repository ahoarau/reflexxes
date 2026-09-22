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

#ifndef INTRINSIC_ICON_REFLEXXES_INTERNAL_MOTION_STATE_WITH_POLYS_H_
#define INTRINSIC_ICON_REFLEXXES_INTERNAL_MOTION_STATE_WITH_POLYS_H_

#include "intrinsic/icon/reflexxes/internal/motion_state.h"
#include "intrinsic/icon/reflexxes/internal/motion_state_transforms.h"
#include "intrinsic/icon/reflexxes/internal/synchronize_dofs_base.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

template <typename T>
struct MotionStateWithPolys;

// This function takes a function and a state object.  The function is executed
// on the state, but in a way that will now generate polynomials to the
// state'state polynomial buffer for every motion segment. Returns the resulting
// application of f(state).
constexpr auto GeneratePolynomials =
    Curry<IsTransform, IsMotionState>([](const auto& f, const auto& state) {
      using S = std::decay_t<decltype(state)>;
      return S(f(MotionStateWithPolys<S>(state)));
    });

// Implementation below

// This is a CRTP wrapper class that allows us to tag a state as "generating
// polynomials".  Any motion application after this is applied will use the
// specialized version below instead of the normal version which will generate
// polynomials.
template <typename T>
struct MotionStateWithPolys : T {
  using T::T;
  explicit MotionStateWithPolys(const T& other) : T(other) {}
};

// This version is specialized to generate polynomials for any type wrapped with
// the MotionStateWithPolys wrapper.
template <typename T>
struct ApplyMotionToState<MotionStateWithPolys<T>> {
  static MotionStateWithPolys<T> Apply(const double target_p,
                                       const double target_v,
                                       const double target_a,
                                       const JerkConstants& j_const,
                                       const double delta_t,
                                       const MotionStateWithPolys<T>& state) {
    double flip_mult = IsFlipped(state) ? -1.0 : 1.0;
    GetPolynomials(state)->TrimSegmentCount(GetPolynomialIndex(state));
    GetPolynomials(state)->AddSegment({.position = flip_mult * GetP(state),
                                       .velocity = flip_mult * GetV(state),
                                       .acceleration = flip_mult * GetA(state),
                                       .jerk = flip_mult * j_const.full,
                                       .start_time = GetT(state),
                                       .end_time = GetT(state) + delta_t});
    return SetPVAAndIncrementT(
        target_p, target_v, target_a, delta_t,
        SetPolynomialIndex(GetPolynomials(state)->GetSegmentCount(), state));
  }
};

// This specialization prevents speculative branches of state from writing out
// polynomials, which is a waste as they are nuked after the speculation.
template <typename T>
struct SetSpeculativeHelper<MotionStateWithPolys<T>> {
  static auto Apply(const MotionStateWithPolys<T>& state) { return T(state); }
};

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_INTERNAL_MOTION_STATE_WITH_POLYS_H_
