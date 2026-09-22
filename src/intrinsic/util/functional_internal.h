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

#ifndef INTRINSIC_UTIL_FUNCTIONAL_INTERNAL_H_
#define INTRINSIC_UTIL_FUNCTIONAL_INTERNAL_H_

#include <cstddef>
#include <functional>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>

namespace intrinsic {
namespace functional {
namespace details {

// The base "arg validator" type.  Acts as a tag for types to indicate they
// validate arguments.
struct ArgValidator {};

// ArgValidator for plain types.  Just asserts that the test type is convertible
// to the type represented by this validator.
template <typename ToType>
struct TypeValidator : ArgValidator {
  template <typename TypeToValidate>
  constexpr static bool Validate() {
    constexpr bool result = std::is_convertible_v<TypeToValidate, ToType>;
    static_assert(result,
                  "ArgValidation failed for curried argument: Argument does "
                  "not convert to the right type.");
    return result;
  }
};

// The "wild card" validator - allows any type to pass.
struct AllowAny : ArgValidator {
  template <typename TypeToValidate>
  constexpr static bool Validate() {
    return true;
  }
};

// This struct and its specialization provides a switch between "plain types"
// and ArgValidator types, such that the user can either pass a plain type and
// have it wrapped by TypeValidator<T> or an ArgValidator which is passed
// through unchanged.
template <typename T, typename Enable = void>
struct ArgValidatorGenerator {
  using type = TypeValidator<T>;
};

template <typename T>
struct ArgValidatorGenerator<
    T, std::enable_if_t<std::is_base_of_v<ArgValidator, T>>> {
  using type = T;
};

// Given a mixed list of types and ArgValidators, produce a tuple of those types
// after being passed through ArgValidatorGenerator to produce a tuple of only
// ArgValidators
template <typename... Ts>
using ArgValidationTuple =
    std::tuple<typename ArgValidatorGenerator<std::decay_t<Ts>>::type...>;

// These templates generate a tuple of "AllowAny" size N, for use in
// ArgValidationTupleN below
template <size_t>
using AnyAtIndex = AllowAny;

// For each index, create a type of AllowAny, thus producing a
// tuple<AllowAny...> size sizeof...(Is)
template <size_t... Is>
auto AnyTuple(std::index_sequence<Is...>) {
  return std::tuple<AnyAtIndex<Is>...>();
}

// Return a tuple of size N with all AllowAny types.
template <size_t N>
auto AnyTuple() {
  return AnyTuple(std::make_index_sequence<N>());
}

// Produce a tuple of size N, with the first portion being arg validators for
// Ts, and the second (N - sizeof...(Ts)) being AllowAny validators.
template <size_t N, typename... Ts>
using ArgValidationTupleN = decltype(std::tuple_cat(
    std::declval<ArgValidationTuple<Ts...>>(), AnyTuple<N - sizeof...(Ts)>()));

// The primary implementation class for the Curry type.  This class keeps
// track of a few things:
//   * The function to be called
//   * A set of arg validator types, one for each argument to be passed to the
//   tuple
//   * The set of bound args, copied as a tuple.
template <typename FunctorT, typename ArgValidationTuple, typename... BoundArgs>
struct CurriedFunctor {
  constexpr static size_t kTotalArgCount =
      std::tuple_size_v<ArgValidationTuple>;

  template <typename F>
  constexpr CurriedFunctor(F&& f, std::tuple<BoundArgs...>&& bound_args)
      : f(std::forward<F>(f)),
        bound_args(std::forward<std::tuple<BoundArgs...>>(bound_args)) {}

  template <size_t Index, typename... Ts>
  constexpr static bool ArgIsValid() {
    return std::tuple_element_t<Index + sizeof...(BoundArgs),
                                ArgValidationTuple>::
        template Validate<std::tuple_element_t<Index, std::tuple<Ts...>>>();
  }

  template <size_t Index, typename... Ts>
  constexpr static bool ValidateArgs() {
    if constexpr (Index < sizeof...(Ts)) {
      constexpr bool is_valid = ArgIsValid<Index, Ts...>();
      static_assert(is_valid,
                    "Arg did not pass evaluation.  See the first parameter of "
                    "this template for index.");
      return is_valid && ValidateArgs<Index + 1, Ts...>();
    }
    return true;
  }

  // The actual work of this class
  // If we have all the args we need (kTotalArgCount) then just call the
  // function.  Otherwise create a new functor with the additional passed in
  // arguments.
  template <typename... Ts>
  constexpr auto operator()(Ts&&... ts) const {
    static_assert(sizeof...(ts) > 0 &&
                      sizeof...(ts) + sizeof...(BoundArgs) <= kTotalArgCount,
                  "Too many arguments to CurriedFunctor operator");
    static_assert(ValidateArgs<0, std::decay_t<Ts>...>());

    if constexpr (sizeof...(Ts) == kTotalArgCount) {
      // All args passed in, this is just a function call on f
      return f(std::forward<Ts>(ts)...);
    } else {
      auto next_args =
          std::tuple_cat(bound_args, std::make_tuple(std::forward<Ts>(ts)...));

      if constexpr (kTotalArgCount > sizeof...(BoundArgs) + sizeof...(ts)) {
        // Generate a new curried functor with the additional args
        return CurriedFunctor<FunctorT, ArgValidationTuple, BoundArgs...,
                              std::decay_t<Ts>...>(f, std::move(next_args));
      } else {
        // Call with our bounded args and the new passed in args
        return std::apply(f, std::move(next_args));
      }
    }
  }

  FunctorT f;
  std::tuple<BoundArgs...> bound_args;
};

// This helper is used to automatically deduce argument types for Curry from
// standard function types - that is function pointers and std::function.  The
// underlying Curry functor will be constructed with the argument parameters set
// to the deduced args of the function.
template <typename F>
struct CurriedFunctorGenerator {};

template <typename F>
using CurriedFunctorT = typename CurriedFunctorGenerator<std::decay_t<F>>::type;

template <typename R, typename... Args>
struct CurriedFunctorGenerator<R (*)(Args...)> {
  using type = CurriedFunctor<R (*)(Args...), ArgValidationTuple<Args...>>;
};

template <typename R, typename... Args>
struct CurriedFunctorGenerator<std::function<R(Args...)>> {
  using type =
      CurriedFunctor<std::function<R(Args...)>, ArgValidationTuple<Args...>>;
};

// Provide an "is_callable" type trait, which returns true if a type is a
// standard "function" type or is a class with an operator() method (any
// operator(), we don't care about which one). This is accomplished by creating
// a type `TestFunctor` that inherits from T and provides an additional
// operator(). If T has any operator() already, this will force an overload. We
// can then test `IsOverloaded` with `TestFunctor`. In the overloaded case,
// `decltype(&TestFunctor::operator())` will fail to resolve when trying to take
// the address of the method.  SFINAE will then route to the "true" path for
// which this decltype does not compile, and thus operator() is overloaded.
// False is where decltype is allowable, and thus the only operator() method is
// the one added by `TestFunctor`.
template <typename T, typename D = std::decay_t<T>, bool = std::is_class_v<D>>
struct is_callable : std::is_function<std::remove_pointer_t<D>> {};

template <typename T>
struct is_callable<T, std::decay_t<T>, true> {
 private:
  struct TestFunctor : std::decay_t<T>, std::function<void()> {};

  template <typename, typename = void>
  struct IsOverloaded : std::true_type {};

  template <typename F>
  struct IsOverloaded<F, std::void_t<decltype(&F::operator())>>
      : std::false_type {};

 public:
  static constexpr bool value = IsOverloaded<TestFunctor>::value;
};

template <typename T>
constexpr bool is_callable_v = is_callable<T>::value;

}  // namespace details
}  // namespace functional
}  // namespace intrinsic

#endif  // INTRINSIC_UTIL_FUNCTIONAL_INTERNAL_H_
