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

// This file provides a variety of helper functionts that assist with
// programming in a functional manner by providing several classic idioms for
// use with C++ functions.
//
// All provided functions are templated, and designed to work with all C++
// callable types.
//
// Compose & Pipe
// One key method of working with functional code is via function composition.
// In this methodolgy, functions are composed into functional chains, where the
// output from one function is passed to the input of the next function in the
// chain.  All functions in the chain (with the exception of the initial call)
// must be unary, and their output must be convertible to the input of the next
// function.
// The two version perform the same operation, where compose calls the functions
// in order from right to left, and pipe proceeds in "reading order" left to
// right.  In addition the >> operator can be used to get a "unix pipe"
// style of programming.
//
// Examples:
//    constexpr auto add_1 = [](int x) { return x + 1; };
//    constexpr auto mult_2 = [](int x) { return x * 2; };
//    constexpr auto mult_3 = [](int x) { return x * 3; };
//    // Each of the below is equivalent to calling mult_2(add_1(mult_3(x)))
//    constexpr auto composed_func = Compose(mult_2, add_1, mult_3);
//    constexpr auto piped_func = Compose(mult_3, add_1, mult_2);
//    constexpr auto piped_func2 = mult_3 >> add_1 >> mult_2;
//    printf("%d == %d == %d\n", compsed_func(2), piped_func(2), piped_func2(2);
//
//    Ouput: 14 == 14 == 14
//
// Currying
// Currying is the technique of translating the evaluation of a function that
// takes multiple arguments into evaluating a sequence of functions, each with a
// single argument.  (This is related but not idential to the concept of partial
// application, which can be achieved with std::bind_front.)  Calling "Curry" on
// a function with N args will return a function which can accept 1-N arguments.
// Passing M arguments to the resulting curried functor where M < N will return
// another function that will accept (N-M) arguments.  Calling a curried functor
// with the number of unbound argumnets will execute the original function with
// the bound and unbound arguments.
// Arguments are copied into the resulting functor, resulting in the storage and
// computation cost of copying a tuple with the bound arguments.
//
// Example:
// int add(int x, int, y, int z) { return x + y + z; }
// Curry(add)(1)(2)(3); // 6
// Curry(add)(1)(2, 3); // 6
// Curry(add)(1, 2)(3); // 6
// Curry(add)(1, 2, 3); // 6
//
// Currying Argument Validation
// The Curry function "lazily" binds the arguments to the call.  The caller is
// not required to pass a specific overload of a function to the Curry function.
// The one requirement it has is that it is provided the number of arguments
// ahead of time (so it can determine when a curry call is made vs the actual
// call).  The caller may still wish to have argument types be validated.  In
// the case of functions, this is done automatically:
//
// int foo(int, int);
// constexpr auto curry_foo = Curry(foo);
// constexpr auto foo_1 = curry_foo(1);  // Returns int(int)
// constexpr auto foo_bad = curry_foo("hello"); // Not an int, won't compile
//
// Because the function is well specified (in this case trivial, but could be
// done via an overload), no additional info is necessary for the Curry function
// to "do the right thing". But what about cases where a parameter is templated?
// In that case, the Curry function provides "ArgValidators" that are run and
// will static assert if a bad parameter is passed to a curry function.  (If the
// arguments are invalid for the call, then this will fail trying to invoke the
// call, much in the same manner as std::bind.  The difference here is that this
// is done on the early "currying" calls as well).
//
// Argument validators are specified at the front of the Curry call, in
// parameter order.  Passing a plain type will assert that the argument must
// pass std::is_convertible.  The validator type AllowAny is a wildcard and will
// not assert anything on the argument.  Additional validators are provided as
// well:
//  - IsBaseOf: Any type that is a base of the given type.
//  - ReturnsResult: Validates the type is callable and returns the given type
//  when called with the specified types.
//  - IsOrReturnsResult: Can take a value, or a function that returns the value
//  when passed the specified types.
//
//  A helper "CurryN" is also provided for a common case where the some or all
//  of the arguments are "AllowAny".  Given the leading template param N, the
//  Curry function will assert any given argument validators, and then AllowAny
//  for the remaining unspecified validators.
//
//  Example:
//  constexpr auto add_3 = Curry(AllowAny, AllowAny, AllowAny,
//     [](auto x, auto y, auto z) { return x + y + z});
//  constexpr auto add_3_b = CurryN<3>([](auto x, auto y, auto z) {
//          return x + y + z});
//  constexpr auto add_3_c = CurryN<3, double>([](double x, auto y, auto z) {
//          return x + y + z});
#ifndef INTRINSIC_UTIL_FUNCTIONAL_H_
#define INTRINSIC_UTIL_FUNCTIONAL_H_

#include <tuple>
#include <type_traits>
#include <utility>

#include "intrinsic/util/functional_internal.h"

namespace intrinsic {
namespace functional {

// Return a function that represents the composition of the given functions,
// such that for f(x)->y and g(y)->z, compose(g,f) returns h(x)->g(f(x))->z.
// All functions but the last one must be unary.
// This call also allows to pass values in as the final argument.  In that case
// the result of the function taking the second parameter as an argument is
// immediately returned.
// The function objects are copied, so keep that in mind if that will incur a
// cost.
template <class F, typename G, typename... Gs>
constexpr auto Compose(F&& f, G&& g, Gs&&... gs) {
  if constexpr (sizeof...(gs) == 0) {
    if constexpr (details::is_callable_v<std::decay_t<G>>) {
      return [f, g](auto&&... xs) {
        return f(g(std::forward<decltype(xs)>(xs)...));
      };
    } else {
      // Just call f with g as the argument.
      return f(g);
    }
  } else {
    return Compose(std::forward<F>(f),
                   Compose(std::forward<G>(g), std::forward<Gs>(gs)...));
  }
}

// Return a function that represents the inverse composition of the given
// functions, such that for f(x)->y and g(y)->z, Pipe(f,g) returns
// h(x)->g(f(x))->z
template <class F, typename G, typename... Gs>
constexpr auto Pipe(F&& f, G&& g, Gs&&... gs) {
  if constexpr (sizeof...(gs) == 0) {
    return Compose(std::forward<G>(g), std::forward<F>(f));
  } else {
    return Pipe(Pipe(std::forward<F>(f), std::forward<G>(g)),
                std::forward<Gs>(gs)...);
  }
}

// TODO(rbutterfoss): Deprecated, remove.
// Convenience operator that pipes two parameters.
template <typename F, typename G>
constexpr auto operator>>(F&& f, G&& g) {
  return Pipe(std::forward<F>(f), std::forward<G>(g));
}

// Convenience operator that pipes two parameters.
template <typename F, typename G>
constexpr auto operator|(F&& f, G&& g) {
  return Pipe(std::forward<F>(f), std::forward<G>(g));
}

// Curry the given function and return its functor.
// If the function is a std::function or function pointer, the arguments are
// deduced automatically.  Otherwise, the template parameters must be given N
// ArgValidators equal to the number of arguments that need to be passed to F.
// If all arguments are to be "unvalidated" and passed AllowAny, the convenience
// function CurryN below may be used.
// See the file header for more info on Argument Validation.
template <typename... ArgValidators, typename F>
constexpr auto Curry(F&& f) {
  if constexpr (sizeof...(ArgValidators) > 0) {
    return details::CurriedFunctor<
        F, details::ArgValidationTuple<ArgValidators...>>(std::forward<F>(f),
                                                          std::tuple<>{});
  } else {
    // Try and get the argument details from the function itself.
    return details::CurriedFunctorT<F>(std::forward<F>(f), std::tuple<>{});
  }
}

// Currys a function F with ArgCount arguments and return its functor.
// If sizeof...(ArgValidators) is < ArgCount, AllowAny will be inserted for the
// remaining unspecified arguments.
// See the file header for more info on Argument Validation.
template <int ArgCount, typename... ArgValidators, typename F>
constexpr auto CurryN(F&& f) {
  return details::CurriedFunctor<
      F, details::ArgValidationTupleN<ArgCount, ArgValidators...>>(
      std::forward<F>(f), std::tuple<>{});
}

// Argument validators for currying

// AllowAny specifies that any argument may work - essentially deferring to the
// compiler to resolve things once the call is made.
using details::AllowAny;

// Validates that the given type is the base of another type.
template <typename Base>
struct IsBaseOf : details::ArgValidator {
  template <typename TypeToValidate>
  constexpr static bool Validate() {
    constexpr bool result = std::is_base_of_v<Base, TypeToValidate>;
    static_assert(result,
                  "ArgValidation failed for curried argument: Argument does "
                  "not have the right base type.");
    return result;
  }
};

// Validates that the given type results in a type convertible to R when called
// with Args
template <typename R, typename... Args>
struct ReturnsResult : details::ArgValidator {
  template <typename TypeToValidate>
  constexpr static bool Validate() {
    constexpr bool result = std::is_invocable_r_v<R, TypeToValidate, Args...>;
    static_assert(result,
                  "ArgValidation failed for curried argument: Argument cannot "
                  "be called properly.");
    return result;
  }
};

// Validates that R is either the result of calling TestType(Args...) or
// TestType itself is convertible to R.
template <typename R, typename... Args>
struct IsOrReturnsResult : details::ArgValidator {
  template <typename TypeToValidate>
  constexpr static bool Validate() {
    if constexpr (std::is_invocable_r_v<R, TypeToValidate, Args...>) {
      return true;
    } else {
      // Default to just seeing if the object type is valid
      return details::TypeValidator<R>::template Validate<TypeToValidate>();
    }
  }
};

// If d, then f(x), else g(x)
constexpr auto IfElse = CurryN<4, bool>(
    [](const bool decision, const auto& f, const auto& g, auto&& x) {
      if (decision) {
        return f(std::forward<decltype(x)>(x));
      } else {
        return g(std::forward<decltype(x)>(x));
      }
    });

}  // namespace functional
}  // namespace intrinsic

#endif  // INTRINSIC_UTIL_FUNCTIONAL_H_
