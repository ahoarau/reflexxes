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

#ifndef INTRINSIC_ICON_REFLEXXES_INTERNAL_POLYNOMIAL_SOLVERS_H_
#define INTRINSIC_ICON_REFLEXXES_INTERNAL_POLYNOMIAL_SOLVERS_H_

#include <array>

#include "intrinsic/util/fixed_vector.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {

// Evaluate an N-th degree polynomial with the given value x.
// Coefficients are specified from least degree to highest degree.
template <size_t N>
double EvalPolynomial(const std::array<double, N>& coefficients,
                      const double x) {
  double result = 0.0;
  double pow = 1.0;
  for (size_t i = 0; i < N; ++i) {
    result += coefficients[i] * pow;
    pow *= x;
  }
  return result;
}

// The biggest number of possible roots from any root solving function is 4 for
// a degree 4 polynomial.
using PolynomialRoots = FixedVector<double, 4>;

// Returns the set of real roots from a degree 1 polynomial.
// Coefficients are specified from least degree to highest degree.
PolynomialRoots CalculatePolynomialRoots(
    const std::array<double, 2>& coefficients);

// Returns the set of real roots from a degree 2 polynomial.
// Coefficients are specified from least degree to highest degree.
PolynomialRoots CalculatePolynomialRoots(
    const std::array<double, 3>& coefficients);

// Always returns a set of 2 real roots by zeroing out the discriminant if it is
// less than zero.  If the coefficients would return 0 or 1 roots, then 0's are
// appended to ensure result is always size 2.
PolynomialRoots SafeCalculatePolynomialRoots(
    const std::array<double, 3>& coefficients);

// Returns the set of real roots from a degree 3 polynomial.
// Coefficients are specified from least degree to highest degree.
PolynomialRoots CalculatePolynomialRoots(
    const std::array<double, 4>& coefficients);

// Returns the set of real roots from a degree 4 polynomial.
// Coefficients are specified from least degree to highest degree.
PolynomialRoots CalculatePolynomialRoots(
    const std::array<double, 5>& coefficients);

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_INTERNAL_POLYNOMIAL_SOLVERS_H_
