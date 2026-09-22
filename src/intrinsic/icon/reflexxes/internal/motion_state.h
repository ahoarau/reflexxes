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

// This file represents the concept of a Motion State.  A Motion State
// represents the current state of the calculations that the decision trees and
// profile methods perform as the central function of the Reflexxes algorithm.
// The state holds not only physical variables like the current acceleration and
// amount of time passed, but also variables like the current profile selected,
// or debugging traits like the current trace of decisions made in a decision
// tree.
// The State is immutable, so functions that deal in state will take a const& to
// the MotionState and return a value to a new MotionState.  Because of this,
// the MotionState objects are constructed in such a way to make the cost of
// doing so in and out of functions to be as inexpensive as possible.  The true
// "state" consists of many many variables but the actual state is represented
// by a small structure that contains just the variables that can mutate, as
// well as a pointer to a StateBase object that contains all the "read only"
// variables that do not change during the course of an evaluation where a
// MotionState is passed around.
// The details about the underlying structure of the state are just that:
// details.  In practice, a user of the state object need only concern
// themselves with the accessor and Mutator functions defined in this file, as
// well as the constructors of the state object.  The only interaction to the
// object should be via these types - the MotionState object itself should be
// thought of as opaque and subject to change.

#ifndef INTRINSIC_ICON_REFLEXXES_INTERNAL_MOTION_STATE_H_
#define INTRINSIC_ICON_REFLEXXES_INTERNAL_MOTION_STATE_H_
#include <cstdint>

#include "intrinsic/icon/reflexxes/constants_external.h"
#include "intrinsic/icon/reflexxes/inputs.h"
#include "intrinsic/icon/reflexxes/internal/constants.h"
#include "intrinsic/icon/reflexxes/motion_polynomials.h"
#include "intrinsic/icon/reflexxes/outputs.h"
#include "intrinsic/icon/reflexxes/profile.h"
#include "intrinsic/util/functional.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

using intrinsic::functional::Curry;

// Constant values for use around a single value of "jerk".
struct JerkConstants {
  explicit JerkConstants(const double j)
      : full(j),
        dbl(j * 2.0),
        half(j * 0.5),
        inv(1.0 / j),
        neg_inv(-inv),
        half_inv(0.5 * inv),
        neg_half_inv(-half_inv),
        one_sixth(j / 6.0) {}

  const double full;
  const double dbl;
  const double half;
  const double inv;
  const double neg_inv;
  const double half_inv;
  const double neg_half_inv;
  const double one_sixth;
};

// This struct holds all the "base" information that a state uses.  These are
// fixed references to memory that does not change, and therefore can be held in
// a single location within the state.  This shrinks the over all size of the
// State variable to a single reference for a whole bunch of fields, which
// improves performance when state is eventually copied.
// The references passed into this struct are not copied and must outlive this
// object.
struct StateBase {
  // These values are fixed for a specific set of input for a specific DOF
  // during an execution of the Online Trajectory Generation algorithm.  They
  // are grouped here for easy lookup from the state variable.
  struct FixedValues {
    explicit FixedValues(const Inputs::DOF& dof)
        : selected(dof.selected),
          ptrgt(dof.target_position),
          vtrgt(dof.target_velocity),
          vmax(dof.max_velocity),
          vmin(dof.min_velocity),
          amax(dof.max_acceleration),
          amin(dof.min_acceleration),
          jmax(dof.max_jerk),
          jmin(dof.min_jerk),
          jmax_consts(jmax),
          jmin_consts(jmin),
          jerk_constant((-jmin) * jmax / ((-jmin) + jmax)),
          inv_jmax_sub_inv_jmin(jmax_consts.inv - jmin_consts.inv),
          inv_jerk_constant(1.0 / jerk_constant),
          dbl_jerk_constant(2 * jerk_constant),
          amax_to_0_dt(amax / -jmin),
          amax_to_0_dv(amax * amax * jmin_consts.neg_half_inv),
          amin_to_0_dt(-amin / jmax),
          amin_to_0_dv(-amin * amin * jmax_consts.half_inv),
          zero_to_amin_dv(amin * amin * jmin_consts.neg_half_inv) {}

    const bool selected;

    const double ptrgt;

    const double vtrgt;
    const double vmax;
    const double vmin;

    const double amax;
    const double amin;

    const double jmax;
    const double jmin;

    const JerkConstants jmax_consts;
    const JerkConstants jmin_consts;
    const double jerk_constant;
    const double inv_jmax_sub_inv_jmin;
    const double inv_jerk_constant;
    const double dbl_jerk_constant;

    const double amax_to_0_dt;
    const double amax_to_0_dv;
    const double amin_to_0_dt;
    const double amin_to_0_dv;
    const double zero_to_amin_dv;
  };

  StateBase(const Inputs::DOF& dof_input, Outputs::DOF& dof_output,
            Outputs::DOF::SubStep& dof_substep,
            MotionPolynomials* polynomials = nullptr, double tsync = 0.,
            bool get_into_boundaries_fast = false)
      : fixed_values(dof_input),
        fixed_values_flipped(FlipInputParameters(dof_input)),
        dof_output(dof_output),
        substep_output(dof_substep),
        polynomials(polynomials),
        tsync(tsync),
        get_into_boundaries_fast(get_into_boundaries_fast) {}

  const FixedValues fixed_values;
  const FixedValues fixed_values_flipped;
  Outputs::DOF& dof_output;
  Outputs::DOF::SubStep& substep_output;
  MotionPolynomials* const polynomials;
  const double tsync;
  const bool get_into_boundaries_fast;
};

namespace details {
// Enums naming the indices of the State tuple.
enum MotionStateMembers {
  kStateBase = 0,
  kCurrentFixedValues = 1,
  kPosition = 2,
  kVelocity = 3,
  kAcceleration = 4,
  kTime = 5,
  kResult = 6,
  kProfile = 7,
  kDecisionTraceIndex = 8,
  kProfileTraceIndex = 9,
  kPolynomialIndex = 10,
  kMaxBaseProperties
};

// The base tuple for MotionState.
// The template arguments allow for additional types to be added to the Tuple.
// The types as listed here MUST be aligned to the enums above.
template <typename... AdditionalValues>
using MotionStateBaseTuple =
    std::tuple<const StateBase&, const StateBase::FixedValues*, const double,
               const double, const double, const double, const bool,
               const Profile, const uint8_t, const uint8_t, const uint8_t,
               AdditionalValues...>;

// Base "tag" struct to compare against.
struct MotionStateBase {};
}  // namespace details

// The base motion state type, which is the above tuple (with any additional
// fields) and the tag type that indicates this type is a motion state.
template <typename... AdditionalValues>
struct MotionStateT : details::MotionStateBaseTuple<AdditionalValues...>,
                      details::MotionStateBase {
  using details::MotionStateBaseTuple<
      AdditionalValues...>::MotionStateBaseTuple;
};

template <typename... AdditionalValues>
MotionStateT<AdditionalValues...> InitializeMotionState(
    const StateBase& state_base, const Inputs::DOF& dof_input,
    AdditionalValues... additional_values) {
  return MotionStateT<AdditionalValues...>(
      state_base, &state_base.fixed_values, dof_input.position,
      dof_input.velocity, dof_input.acceleration, 0, false, Profile::kUndefined,
      0, 0, 0, additional_values...);
}

// A Curry argument validator for "being a motion state"
using IsMotionState = intrinsic::functional::IsBaseOf<details::MotionStateBase>;

// Getters & Setters

// Define a getter for each member of the State tuple.

constexpr auto GetP = [](const auto& state) {
  return std::get<details::MotionStateMembers::kPosition>(state);
};

constexpr auto GetV = [](const auto& state) {
  return std::get<details::MotionStateMembers::kVelocity>(state);
};

constexpr auto GetA = [](const auto& state) {
  return std::get<details::MotionStateMembers::kAcceleration>(state);
};

constexpr auto GetT = [](const auto& state) {
  return std::get<details::MotionStateMembers::kTime>(state);
};

constexpr auto GetResult = [](const auto& state) {
  return std::get<details::MotionStateMembers::kResult>(state);
};

constexpr auto GetProfile = [](const auto& state) {
  return std::get<details::MotionStateMembers::kProfile>(state);
};

constexpr auto GetDecisionTraceIndex = [](const auto& state) {
  return std::get<details::MotionStateMembers::kDecisionTraceIndex>(state);
};

constexpr auto GetDecisionProfileIndex = [](const auto& state) {
  return std::get<details::MotionStateMembers::kProfileTraceIndex>(state);
};

constexpr auto GetPolynomialIndex = [](const auto& state) {
  return std::get<details::MotionStateMembers::kPolynomialIndex>(state);
};

constexpr auto GetFixedValues = [](const auto& state) {
  return std::get<details::MotionStateMembers::kCurrentFixedValues>(state);
};

constexpr auto GetStateBase = [](const auto& state) -> const StateBase& {
  return std::get<details::MotionStateMembers::kStateBase>(state);
};

// Define a getter for each member of the StateBase object under the state
// tuple.

constexpr auto GetDOFOutput = [](const auto& state) -> Outputs::DOF& {
  return GetStateBase(state).dof_output;
};

constexpr auto GetDecisionTrace =
    [](const auto& state) -> Outputs::DOF::SubStep::DecisionTrace& {
  return GetStateBase(state).substep_output.decision_trace;
};

constexpr auto GetProfileTrace =
    [](const auto& state) -> Outputs::DOF::SubStep::ProfileTrace& {
  return GetStateBase(state).substep_output.applied_profile_trace;
};

constexpr auto GetPolynomials = [](const auto& state) {
  return GetStateBase(state).polynomials;
};

constexpr auto GetTSync = [](const auto& state) {
  return GetStateBase(state).tsync;
};

constexpr auto GetIntoBoundariesFast = [](const auto& state) {
  return GetStateBase(state).get_into_boundaries_fast;
};

// Define a getter for each member under the FixedValues struct.
// This is the "current" fixed values - we don't care about the other one under
// StateBase.
constexpr auto GetPTrgt = [](const auto& state) {
  return GetFixedValues(state)->ptrgt;
};

constexpr auto GetVTrgt = [](const auto& state) {
  return GetFixedValues(state)->vtrgt;
};

constexpr auto GetVMax = [](const auto& state) {
  return GetFixedValues(state)->vmax;
};

constexpr auto GetVMin = [](const auto& state) {
  return GetFixedValues(state)->vmin;
};

constexpr auto GetAMax = [](const auto& state) {
  return GetFixedValues(state)->amax;
};

constexpr auto GetAMin = [](const auto& state) {
  return GetFixedValues(state)->amin;
};

constexpr auto GetJMaxConsts = [](const auto& state) {
  return GetFixedValues(state)->jmax_consts;
};

constexpr auto GetJMinConsts = [](const auto& state) {
  return GetFixedValues(state)->jmin_consts;
};

constexpr auto GetJerkConstant = [](const auto& state) {
  return GetFixedValues(state)->jerk_constant;
};

constexpr auto GetInvJerkConstant = [](const auto& state) {
  return GetFixedValues(state)->inv_jerk_constant;
};

constexpr auto GetInvJMaxSubInvJMin = [](const auto& state) {
  return GetFixedValues(state)->inv_jmax_sub_inv_jmin;
};

constexpr auto GetAMaxToZeroDeltaT = [](const auto& state) {
  return GetFixedValues(state)->amax_to_0_dt;
};

constexpr auto GetAMaxToZeroDeltaV = [](const auto& state) {
  return GetFixedValues(state)->amax_to_0_dv;
};

constexpr auto GetAMinToZeroDeltaT = [](const auto& state) {
  return GetFixedValues(state)->amin_to_0_dt;
};

constexpr auto GetAMinToZeroDeltaV = [](const auto& state) {
  return GetFixedValues(state)->amin_to_0_dv;
};

constexpr auto GetZeroToAMinDeltaV = [](const auto& state) {
  return GetFixedValues(state)->zero_to_amin_dv;
};

// Define a getter for each value under jmax/jmin const structures.

const auto GetJMax = [](const auto& state) {
  return GetJMaxConsts(state).full;
};

const auto GetJMin = [](const auto& state) {
  return GetJMinConsts(state).full;
};

const auto GetInvJMax = [](const auto& state) {
  return GetJMaxConsts(state).inv;
};

const auto GetInvJMin = [](const auto& state) {
  return GetJMinConsts(state).inv;
};

const auto GetNegInvJMin = [](const auto& state) {
  return GetJMinConsts(state).neg_inv;
};

const auto GetHalfInvJMax = [](const auto& state) {
  return GetJMaxConsts(state).half_inv;
};

const auto GetNegHalfInvJMin = [](const auto& state) {
  return GetJMinConsts(state).neg_half_inv;
};

const auto GetDblJMax = [](const auto& state) {
  return GetJMaxConsts(state).dbl;
};

const auto GetDblJMin = [](const auto& state) {
  return GetJMinConsts(state).dbl;
};

const auto GetHalfJMax = [](const auto& state) {
  return GetJMaxConsts(state).half;
};

const auto GetHalfJMin = [](const auto& state) {
  return GetJMinConsts(state).half;
};

const auto GetOneSixthJMax = [](const auto& state) {
  return GetJMaxConsts(state).one_sixth;
};

const auto GetOneSixthJMin = [](const auto& state) {
  return GetJMinConsts(state).one_sixth;
};

namespace details {
// Creates a copy of the given tuple and returns a type whose members are
// non-const versions of the original tuples members.
template <typename... Ts>
auto MakeNonConstTuple(const std::tuple<Ts...>& tuple) {
  return std::tuple<std::remove_const_t<Ts>...>(tuple);
}

// Returns a new item with the Nth value set to val.
template <size_t N, typename StateType, typename T>
StateType SetValue(T&& val, const StateType& state) {
  auto temp = MakeNonConstTuple(state);
  std::get<N>(temp) = std::forward<T>(val);
  return StateType(temp);
}

// Helper function that defines a curried mutator on one of the properties of
// the state tuple.
template <typename T, MotionStateMembers PropertyIndex>
constexpr auto MakeStateMutator() {
  return Curry<T, IsMotionState>([](const T new_value, const auto& state) {
    return details::SetValue<PropertyIndex>(new_value, state);
  });
}

}  // namespace details

// Expose a mutator for each mutable property of the state tuple.

constexpr auto SetP = details::MakeStateMutator<double, details::kPosition>();
constexpr auto SetV = details::MakeStateMutator<double, details::kVelocity>();
constexpr auto SetA =
    details::MakeStateMutator<double, details::kAcceleration>();
constexpr auto SetT = details::MakeStateMutator<double, details::kTime>();
constexpr auto SetResult = details::MakeStateMutator<bool, details::kResult>();
constexpr auto SetDecisionTraceIndex =
    details::MakeStateMutator<uint8_t, details::kDecisionTraceIndex>();
constexpr auto SetProfileTraceIndex =
    details::MakeStateMutator<uint8_t, details::kProfileTraceIndex>();
constexpr auto SetPolynomialIndex =
    details::MakeStateMutator<uint8_t, details::kPolynomialIndex>();

// This function is a concession to the debug builds which run slower than the
// reflexxes predecessor because of the functional patterns.  This function is
// the heart of many of the motion transforms (via ApplyMotionToState), so
// allowing these frequent operations to be "pre-optimized" by not using
// functional composition with the setters/getters helps the debug build quite a
// bit.
template <typename StateType>
StateType SetPVAAndIncrementT(const double p, const double v, const double a,
                              const double dt, const StateType& state) {
  auto temp = details::MakeNonConstTuple(state);
  std::get<details::MotionStateMembers::kPosition>(temp) = p;
  std::get<details::MotionStateMembers::kVelocity>(temp) = v;
  std::get<details::MotionStateMembers::kAcceleration>(temp) = a;
  std::get<details::MotionStateMembers::kTime>(temp) += dt;
  return StateType(temp);
}

// Increments t by dt
constexpr auto IncrementT =
    Curry<double, IsMotionState>([](const double dt, const auto& state) {
      return SetT(GetT(state) + dt, state);
    });

// Whether or not we're using the flipped input or not
template <typename State>
bool IsFlipped(const State& state) {
  return GetFixedValues(state) ==
         &std::get<details::MotionStateMembers::kStateBase>(state)
              .fixed_values_flipped;
}

// Sets/Gets whether or not this state is the result of a successful operation.
constexpr auto SetSuccess = [](const auto& state) {
  return SetResult(true, state);
};
constexpr auto IsSuccess = [](const auto& state) { return GetResult(state); };

// Sets/Gets whether or not this state is the result of a failed operation.
constexpr auto SetFailure = [](const auto& state) {
  return SetResult(false, state);
};
constexpr auto IsFailure = [](const auto& state) { return !GetResult(state); };

// Appends a decision (the step and the true/false value of the decision) to the
// given state'state decision trace buffer.
constexpr auto AppendDecisionTrace = Curry<int, bool, IsMotionState>(
    [](int step_num, bool decision, const auto& state) {
      if constexpr (kEnableTracing) {
        int index = GetDecisionTraceIndex(state);
        CHECK_LT(index, kTraceSize);
        GetDecisionTrace(state).AddDecisionTrace(
            index, decision, static_cast<uint8_t>(step_num));
        return SetDecisionTraceIndex(index + 1, state);
      } else {
        return state;
      }
    });

// Appends the given profile to the given state'state profile buffer.
constexpr auto AppendProfileTrace =
    Curry<Profile, IsMotionState>([](const Profile profile, const auto& state) {
      if constexpr (kEnableTracing) {
        int index = GetDecisionProfileIndex(state);
        CHECK_LT(index, kTraceSize);
        GetProfileTrace(state)[index] = profile;
        return SetProfileTraceIndex(index + 1, state);
      } else {
        return state;
      }
    });

// Sets the applied profile of the state, while also appending to the profile
// trace and setting this to "success".
constexpr auto SetProfile =
    Curry<Profile, IsMotionState>([](const Profile profile, const auto& state) {
      return details::SetValue<details::MotionStateMembers::kProfile>(profile,
                                                                      state) |
             AppendProfileTrace(profile) | SetSuccess;
    });

// Flips the input and current motion params
constexpr auto FlipState = [](const auto& state) {
  using intrinsic::functional::operator|;

  const StateBase& base =
      std::get<details::MotionStateMembers::kStateBase>(state);
  return details::SetValue<details::MotionStateMembers::kCurrentFixedValues>(
             IsFlipped(state) ? &base.fixed_values : &base.fixed_values_flipped,
             state) |
         SetV(-GetV(state)) | SetA(-GetA(state)) | SetP(-GetP(state)) |
         AppendProfileTrace(Profile::kFlip);
};

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_INTERNAL_MOTION_STATE_H_
