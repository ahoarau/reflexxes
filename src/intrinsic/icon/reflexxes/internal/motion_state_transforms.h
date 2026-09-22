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

#ifndef INTRINSIC_ICON_REFLEXXES_INTERNAL_MOTION_STATE_TRANSFORMS_H_
#define INTRINSIC_ICON_REFLEXXES_INTERNAL_MOTION_STATE_TRANSFORMS_H_

#include <algorithm>
#include <cstddef>
#include <type_traits>

#include "intrinsic/icon/reflexxes/inputs.h"
#include "intrinsic/icon/reflexxes/internal/math.h"
#include "intrinsic/icon/reflexxes/internal/motion_state.h"
#include "intrinsic/util/functional.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

using intrinsic::functional::operator|;
using intrinsic::functional::AllowAny;
using intrinsic::functional::IfElse;

// Have some validators for ourselves
using DoubleGetter =
    intrinsic::functional::ReturnsResult<double, MotionStateT<>>;
using DoubleOrGetter =
    intrinsic::functional::IsOrReturnsResult<double, MotionStateT<>>;
using IsTransform =
    intrinsic::functional::ReturnsResult<details::MotionStateBase,
                                         MotionStateT<>>;

// The default reaction to set_speculative being called on a type. (See below).
// The default action is a no-op, in that the speculative action runs no
// different than the non-speculative action.
template <typename State>
struct SetSpeculativeHelper {
  static State Apply(const State& state) { return state; }
};

// This function allows the state to be modified in such a way that certain
// calculations can be ignored.  A speculative calculation is one done to get
// some value used to calculate a motion function before the final state is
// calculated.  The "and_then" functions below make use of this, in that they
// calculate the effects of the actions that come after the initial action to
// take into account their additional velocity or acceleration, and then
// performs its own calculation with that in mind.
// By wrapping a function as "speculative", the state can be modified such that
// actions that are only "non-speculative" will take effect.  This can be done
// on a per type specialization using the above Helper functor.
constexpr auto set_speculative = [](const auto& state) {
  return SetSpeculativeHelper<std::decay_t<decltype(state)>>::Apply(state);
};

// This is the core function that applies motion to a state.  It is templated
// for types of state that wish to use the solver framework below to do
// additional things. (e.g. Write polynomials).
//
// This function is in a struct to allow the function to be partially
// specialized for specific types of state.  (see motion_state_with_polys.h).
template <typename State>
struct ApplyMotionToState {
  // Apply the motion variables resulting from the application of the
  // given jerk for dt to state state.
  inline static State Apply(const double target_p, const double target_v,
                            const double target_a,
                            const JerkConstants& /*unused*/, const double dt,
                            const State& state) {
    return SetPVAAndIncrementT(target_p, target_v, target_a, dt, state);
  }
};

// These helpers work with ApplyMotionToState above to apply jerk with various
// known and unknown values about the jerk.  These variations allow operations
// which may already know the outcomes of certain state variables to be more
// efficient.

// Apply jerk until the target_a is met, with a known ending target_v and dt for
// State state.
template <typename State>
inline State JerkVAT(const JerkConstants& j_const, const double target_v,
                     const double target_a, const double dt,
                     const State& state) {
  double p_prime = GetP(state) + GetV(state) * dt +
                   0.5 * GetA(state) * Power2(dt) +
                   j_const.one_sixth * Power3(dt);
  return ApplyMotionToState<State>::Apply(p_prime, target_v, target_a, j_const,
                                          dt, state);
}

// Apply jerk until the target_a is met, with a known ending target_v for the
// operation.
template <typename State>
inline State JerkVA(const JerkConstants& j_const, const double target_v,
                    const double target_a, const State& state) {
  return JerkVAT(j_const, target_v, target_a,
                 (target_a - GetA(state)) * j_const.inv, state);
}

// Apply jerk until the target_a is met, with a known dt for the operation.
template <typename State>
inline State JerkAT(const JerkConstants& j_const, const double target_a,
                    const double dt, const State& state) {
  return JerkVAT(
      j_const,
      GetV(state) + (Power2(target_a) - Power2(GetA(state))) * j_const.half_inv,
      target_a, dt, state);
}

// Apply jerk until the target_a is met.
template <typename State>
inline State JerkA(const JerkConstants& j_const, const double target_a,
                   const State& state) {
  return JerkAT(j_const, target_a, (target_a - GetA(state)) * j_const.inv,
                state);
}

// Take a value f, where F is either an instance of T or a function of type
// F(S)->T.  Return either the value or the result of calling f(state).
template <typename T, typename F, typename S>
T CallOrGet(const F& f, const S& state) {
  if constexpr (std::is_convertible_v<F, T>) {
    return f;
  } else {
    return f(state);
  }
}

// Always return the constant v.
// Can be used as a constant fill in for any place that needs to take a state
// accessor function.
constexpr auto Always = Curry<AllowAny, IsMotionState>(
    [](const auto value, const auto&) { return value; });

// A no-op method that just returns itself
constexpr auto Identity = [](const auto& value) { return value; };

// Acceleration to the value given by a_or_func/a_or_func(state). (a->a' PosLin)
constexpr auto AUpToAPrime = Curry<DoubleOrGetter, IsMotionState>(
    [](const auto& a_or_func, const auto& state) {
      return JerkA(GetJMaxConsts(state), CallOrGet<double>(a_or_func, state),
                   state);
    });

// Acceleration to the value given by a_or_func/a_or_func(state). (a->a' Neglin)
constexpr auto ADownToAPrime = Curry<DoubleOrGetter, IsMotionState>(
    [](const auto& a_or_func, const auto& state) {
      return JerkA(GetJMinConsts(state), CallOrGet<double>(a_or_func, state),
                   state);
    });

// More acceleration functions to 0, min, max

constexpr auto AUpToAMax = AUpToAPrime(GetAMax);
constexpr auto ADownToAMax = ADownToAPrime(GetAMax);
constexpr auto ADownToAMin = ADownToAPrime(GetAMin);
constexpr auto AUpToZero = AUpToAPrime(0);
constexpr auto ADownToZero = ADownToAPrime(0);

// Since jerk, amax & amin are fixed values, we can precompute their effect on v
// & t when going from amin to 0 and amax to 0.
constexpr auto AMaxDownToZero = [](const auto& state) {
  return JerkVAT(GetJMinConsts(state), GetV(state) + GetAMaxToZeroDeltaV(state),
                 0., GetAMaxToZeroDeltaT(state), state);
};

constexpr auto AMinDownToZero = [](const auto& state) {
  return JerkVAT(GetJMaxConsts(state), GetV(state) + GetAMinToZeroDeltaV(state),
                 0., GetAMinToZeroDeltaT(state), state);
};

// Accelerate to amax and then back to 0. (a.k.a a->amax postri)
constexpr auto AUpToAMaxDownToZero = AUpToAMax | AMaxDownToZero;

// Accelerate to amin and then back to 0. (a.k.a a->amin negtri)
constexpr auto ADownToAMinUpToZero = ADownToAMin | AMinDownToZero;

// Accelerate up to the target_a and then back to 0.
// target_a can be a value or a function double(State).
constexpr auto AUpToAPrimePosTri = Curry<DoubleOrGetter, IsMotionState>(
    [](const auto& a_or_target_a_func, const auto& state) {
      return AUpToAPrime(a_or_target_a_func, state) | ADownToZero;
    });

// Accelerate down to the target a and then back to 0
// target_a can be a value or a function double(State).
constexpr auto ADownToAPrimeNegTri = Curry<DoubleOrGetter, IsMotionState>(
    [](const auto& target_a_or_func, const auto& state) {
      return ADownToAPrime(target_a_or_func, state) | AUpToZero;
    });

// Hold at the current v, a for the amount of time given by dt_or_func, which
// can be the value or a function double(State) that returns the value.
constexpr auto HoldForDeltaT = Curry<DoubleOrGetter, IsMotionState>(
    [](const auto& dt_or_func, const auto& state) {
      static const JerkConstants zero_jerk(0.);

      double dt = CallOrGet<double>(dt_or_func, state);
      return JerkVAT(zero_jerk, GetV(state) + dt * GetA(state), GetA(state), dt,
                     state);
    });

// Hold until the time given by applying target_t_func to state, accounting for
// the time taken to accomplish "after_func".  after_func is then applied to the
// result of the hold.
constexpr auto HoldToTPrimeAndThen =
    Curry<DoubleGetter, IsTransform, IsMotionState>(
        [](const auto& target_t_func, const auto& after_func,
           const auto& state) {
          return HoldForDeltaT(target_t_func(state) -
                                   GetT(after_func(set_speculative(state))),
                               state) |
                 after_func;
        });

// Hold at current v,a until t=tsync, also accounting for the time in the after
// func, and then perform after_func.
constexpr auto HoldToTSyncAndThen = HoldToTPrimeAndThen(GetTSync);

// Hold at current v,a until t=tsync for a given state state
constexpr auto HoldToTSync = HoldToTPrimeAndThen(GetTSync, Identity);

// Hold until the given change in velocity is met at the current a.
constexpr auto HoldToDeltaV =
    Curry<double, IsMotionState>([](const double dv, const auto& state) {
      return HoldForDeltaT(dv / GetA(state), state);
    });

// Holds until v is vmin
constexpr auto HoldToVMin = [](const auto& state) {
  return HoldToDeltaV(GetVMin(state) - GetV(state), state);
};

// Holds until v->v', while accounting for additional changes to v due to
// after_func. Does not actually apply after_func, just adjusts hold time to
// deal with it.
constexpr auto HoldToVPrimeAccountingFor =
    Curry<DoubleGetter, IsTransform, IsMotionState>(
        [](const auto& target_v_func, const auto& after_func,
           const auto& state) {
          return HoldToDeltaV(
              target_v_func(state) - GetV(after_func(set_speculative(state))),
              state);
        });
constexpr auto HoldToVMaxAccountingFor = HoldToVPrimeAccountingFor(GetVMax);

// Holds until v->v', while accounting for additional changes to v due to
// after_func, and then applies after_func.
constexpr auto HoldToVPrimeAndThen =
    Curry<DoubleGetter, IsTransform, IsMotionState>(
        [](const auto& target_v_func, const auto& after_func,
           const auto& state) {
          return HoldToVPrimeAccountingFor(target_v_func, after_func, state) |
                 after_func;
        });
constexpr auto HoldToVTrgtAndThen = HoldToVPrimeAndThen(GetVTrgt);

// a->+ahld->hold->func such that v=target_v_func(state) and t =
// target_t_func(state)
constexpr auto ADownToAHldThenHoldAndThen =
    Curry<DoubleGetter, DoubleGetter, IsTransform, IsMotionState>(
        [](const auto& target_v_func, const auto& target_t_func,
           const auto& after_func, const auto& state) {
          double t_prime = target_t_func(state);
          auto s_after = set_speculative(state) | after_func;
          if (GetT(s_after) < t_prime) {
            return state |
                   ADownToAPrime((target_v_func(state) - GetV(s_after)) /
                                 (t_prime - GetT(s_after))) |
                   HoldToTPrimeAndThen(target_t_func, after_func);
          }
          return after_func(state);
        });

// Brings a ->ahld->hold->after_func such that v=vtrgt and t=tsync
constexpr auto ADownToAHldThenHoldToVTrgtTsyncAndThen =
    ADownToAHldThenHoldAndThen(GetVTrgt, GetTSync);

// v->v' via Postrap accounting for after_func, and then applying after_func.
// Assumes v < v`
constexpr auto VUpToVPrimePosTrapAndThen =
    Curry<DoubleGetter, IsTransform, IsMotionState>(
        [](const auto& target_v_func, const auto& after_func,
           const auto& state) {
          return AUpToAMax(state) |
                 HoldToVPrimeAndThen(target_v_func,
                                     AMaxDownToZero | after_func);
        });
constexpr auto VUpToVTrgtPosTrapAndThen = VUpToVPrimePosTrapAndThen(GetVTrgt);
constexpr auto VUpToVTrgtPosTrap = VUpToVTrgtPosTrapAndThen(Identity);
constexpr auto VUpToVMaxPosTrap = VUpToVPrimePosTrapAndThen(GetVMax, Identity);
constexpr auto VUpToVMinPosTrap = VUpToVPrimePosTrapAndThen(GetVMin, Identity);

// v->v' via Negtrap(a->amin->hold->0)
// Assumes v > v`
constexpr auto VDownToVPrimeNegTrap = Curry<DoubleGetter, IsMotionState>(
    [](const auto& target_v_func, const auto& state) {
      return ADownToAMin(state) |
             HoldToVPrimeAndThen(target_v_func, AMinDownToZero);
    });
constexpr auto VDownToVTrgtNegTrap = VDownToVPrimeNegTrap(GetVTrgt);
constexpr auto VDownToVMinNegTrap = VDownToVPrimeNegTrap(GetVMin);

// v->v' by a (up)->a' such that v=v' at a=a'.
// Assumes v < v`
constexpr auto VUpToVPrimePosLin = Curry<DoubleGetter, IsMotionState>(
    [](const auto target_v_func, const auto& state) {
      auto target_a =
          GetSqrt(Power2(GetA(state)) +
                  GetDblJMax(state) * (target_v_func(state) - GetV(state)));
      return AUpToAPrime(target_a, state);
    });
constexpr auto VUpToVMaxPosLin = VUpToVPrimePosLin(GetVMax);
constexpr auto VUpToVMinPosLin = VUpToVPrimePosLin(GetVMin);

// v->v' by a(down)->a' such that v=v' at a=a'.
constexpr auto VDownToVPrimeNegLin = Curry<DoubleGetter, IsMotionState>(
    [](const auto& target_v_func, const auto& state) {
      auto target_a =
          GetSqrt(Power2(GetA(state)) +
                  GetDblJMin(state) * (target_v_func(state) - GetV(state)));
      return ADownToAPrime(target_a, state);
    });
constexpr auto VDownToVMinNegLin = VDownToVPrimeNegLin(GetVMin);

// Calculate the peak acceleration (apeak) to meet the velocity given by
// target_v or target_v(state) after accelerating a(up)->apeak->0
constexpr auto CalcPeakAccelVPrimePosTri = Curry<DoubleOrGetter, IsMotionState>(
    [](const auto& target_v_or_func, const auto& state) {
      return GetSqrt(
          ((Power2(GetA(state)) * GetInvJMax(state)) +
           2.0 * (CallOrGet<double>(target_v_or_func, state) - GetV(state))) *
          GetJerkConstant(state));
    });

// Calculate apeak (PosTri) to hit vmin for a given state.
constexpr auto CalcPeakAccelVMinPosTri = CalcPeakAccelVPrimePosTri(GetVMin);

// Calculate apeak (PosTri) to hit vtrgt for a given state.
constexpr auto CalcPeakAccelVTrgtPosTri = CalcPeakAccelVPrimePosTri(GetVTrgt);

// Calculate apeak (PosTri) to hit vmax for a given state.
constexpr auto CalcPeakAccelVMaxPosTri = CalcPeakAccelVPrimePosTri(GetVMax);

// Accelerate PosTri such that v is vmin afterwards.
// Assumes v < vmin
constexpr auto VUpToVMinPosTri = AUpToAPrimePosTri(CalcPeakAccelVMinPosTri);

// Accelerate PosTri such that v is vtrgt afterwards.
// Assumes v < vtrgt
constexpr auto VUpToVTrgtPosTri = AUpToAPrimePosTri(CalcPeakAccelVTrgtPosTri);

// Accelerate PosTri such that v is vmax afterwards.
// Assumes v < vmax
constexpr auto VUpToVMaxPosTri = AUpToAPrimePosTri(CalcPeakAccelVMaxPosTri);

// Perform a->apeak->0 and then after_func with apeak such that v =
// target_vf(state) afterwards Assumes v < vmin
constexpr auto VUpToVPrimePosTriAndThen =
    Curry<DoubleGetter, IsTransform, IsMotionState>(
        [](const auto& target_vf, const auto& after_func, const auto& state) {
          // Calc the dv from the after func, with the knowledge we will hit a=0
          double after_dv =
              GetV(after_func(set_speculative(state) | SetA(0))) - GetV(state);
          return AUpToAPrimePosTri(
                     CalcPeakAccelVPrimePosTri(target_vf(state) - after_dv),
                     state) |
                 after_func;
        });
constexpr auto VToVTrgtPosTriAndThen = VUpToVPrimePosTriAndThen(GetVTrgt);

// Calculate the peak acceleration (apeak) to meet the velocity given by
// target_v or target_v(state) after accelerating a(down)->apeak->0
constexpr auto CalcPeakAccelVPrimeNegTri = Curry<DoubleGetter, IsMotionState>(
    [](const auto& target_v_func, const auto& state) -> double {
      return -1.0 * GetSqrt(((Power2(GetA(state)) / (-GetJMin(state))) -
                             2.0 * (target_v_func(state) - GetV(state))) *
                            GetJerkConstant(state));
    });
constexpr auto VDownToVTrgtNegTri =
    ADownToAPrimeNegTri(CalcPeakAccelVPrimeNegTri(GetVTrgt));
constexpr auto VDownToZeroNegTri =
    ADownToAPrimeNegTri(CalcPeakAccelVPrimeNegTri(Always(0)));

// Holds until p given by target_pf is met, assuming the difference in change
// in p brought by after_func, which is applied after the hold.
constexpr auto HoldToPPrimeAndThen =
    Curry<DoubleGetter, IsTransform, IsMotionState>(
        [](const auto& target_pf, const auto& after_func, const auto& state) {
          return HoldForDeltaT((target_pf(state) -
                                GetP(after_func(set_speculative(state)))) /
                                   GetV(state),
                               state) |
                 after_func;
        });
constexpr auto HoldToPTrgtAndThen = HoldToPPrimeAndThen(GetPTrgt);

// Ensures v is no larger than vmax
constexpr auto CapVToVMax = [](const auto& state) {
  return SetV(std::min(GetV(state), GetVMax(state)), state);
};

// State Comparisons
// These comparison functions are used inline for functions that return a bool.
// For decision trees, a common action is to calculate a "what-if" and evaluate
// some feature of the state and see how it relates to some target feature (ex:
// "if acceleration is set to max and then 0 at max jerk, is position less than
// the target position?").  The use of these functions allow a pure functional
// way to write the decisions.
// The naming scheme is Cmp(Attribute)(Op)(TargetAttribute), so for the above
// you would use CmpPLTPTrgt to compare the two.

constexpr auto CmpGt = Curry<DoubleGetter, DoubleGetter, IsMotionState>(
    [](const auto& fa, const auto& fb, const auto& state) {
      return fa(state) > fb(state);
    });
constexpr auto CmpGte = Curry<DoubleGetter, DoubleGetter, IsMotionState>(
    [](const auto& fa, const auto& fb, const auto& state) {
      return fa(state) >= fb(state);
    });
constexpr auto CmpLt = Curry<DoubleGetter, DoubleGetter, IsMotionState>(
    [](const auto& fa, const auto& fb, const auto& state) {
      return fa(state) < fb(state);
    });
constexpr auto CmpLte = Curry<DoubleGetter, DoubleGetter, IsMotionState>(
    [](const auto& fa, const auto& fb, const auto& state) {
      return fa(state) <= fb(state);
    });
constexpr auto CmpSignEq = Curry<DoubleGetter, DoubleGetter, IsMotionState>(
    [](const auto& fa, const auto& fb, const auto& state) {
      return GetSign(fa(state)) == GetSign(fb(state));
    });

constexpr auto CmpAGteZero = CmpGte(GetA, Always(0));
constexpr auto CmpALteAMax = CmpLte(GetA, GetAMax);

constexpr auto CmpVGteZero = CmpGte(GetV, Always(0));

constexpr auto CmpVLtVMin = CmpLt(GetV, GetVMin);
constexpr auto CmpVLteVMin = CmpLte(GetV, GetVMin);

constexpr auto CmpVGtVMin = CmpGt(GetV, GetVMin);

constexpr auto CmpVGteVMin = CmpGte(GetV, GetVMin);

constexpr auto CmpVLtVMax = CmpLt(GetV, GetVMax);

constexpr auto CmpVLteVMax = CmpLte(GetV, GetVMax);

constexpr auto CmpVGtVMax = CmpGt(GetV, GetVMax);

constexpr auto CmpVGteVMax = CmpGte(GetV, GetVMax);

constexpr auto CmpVGtVTrgt = CmpGt(GetV, GetVTrgt);

constexpr auto CmpVGteVTrgt = CmpGte(GetV, GetVTrgt);

constexpr auto CmpVLtVTrgt = CmpLt(GetV, GetVTrgt);

constexpr auto CmpVLteVTrgt = CmpLte(GetV, GetVTrgt);

constexpr auto CmpVSignEqVTrgt = CmpSignEq(GetV, GetVTrgt);

constexpr auto CmpVTrgtLtZero = CmpLt(GetVTrgt, Always(0));

constexpr auto CmpVTrgtLteZero = CmpLte(GetVTrgt, Always(0));

// Position comparators
constexpr auto CmpPGtPTrgt = CmpGt(GetP, GetPTrgt);

constexpr auto CmpPGtePTrgt = CmpGte(GetP, GetPTrgt);

constexpr auto CmpPLtPTrgt = CmpLt(GetP, GetPTrgt);

constexpr auto CmpPLtePTrgt = CmpLte(GetP, GetPTrgt);

constexpr auto ComputeEpsilonPTrgt = [](const auto& state) {
  return GetPTrgt(state) + kAbsoluteStep1B1Epsilon +
         kRelativeStep1B1Epsilon * fabs(GetP(state) - GetPTrgt(state));
};

constexpr auto CmpPLteEpsilonPTrgt = CmpLte(GetP, ComputeEpsilonPTrgt);

constexpr auto CmpPGteEpsilonPTrgt = CmpGte(GetP, ComputeEpsilonPTrgt);

// This check subtracts an epsilon based on the current velocity for position
// comparisons
constexpr auto CmpPLtePTrgtVMaxEpsilon = [](const auto& state) {
  if (GetV(state) > 0.0) {
    return CmpPLtePTrgt(state |
                        SetP(GetP(state) - GetV(state) * kStep2VMaxEpsilon));
  }
  return CmpPLtePTrgt(state);
};

constexpr auto CmpPGtePTrgtVMinEpsilon = [](const auto& state) {
  if (GetV(state) < 0.0) {
    return CmpPGtePTrgt(state |
                        SetP(GetP(state) - GetV(state) * kStep2VMaxEpsilon));
  }
  return CmpPGtePTrgt(state);
};

constexpr auto CmpTGteTSync = CmpGte(GetT, GetTSync);

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_INTERNAL_MOTION_STATE_TRANSFORMS_H_
