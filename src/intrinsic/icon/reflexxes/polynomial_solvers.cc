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

#include "intrinsic/icon/reflexxes/internal/polynomial_solvers.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <utility>

#include "intrinsic/icon/reflexxes/internal/functional.h"
#include "intrinsic/icon/reflexxes/internal/math.h"

namespace intrinsic {
namespace reflexxes {
namespace internal {
namespace {
/// The constant used to check whether a quartic equation is biquadratic.
constexpr double kBiQuadraticScaleFactor = 1.0e-20;

/// The scale factor used on std::numeric_limits<double>::epsilon()
/// to check for duplicate real roots in cubic equation.
constexpr double kCubicDoubleRootsEpsilonScaleFactor = 10.0;

/// The absolute value below which the value is considered as zero.
constexpr double kZeroThreshold = 1.0e-6;

// Solves the quadratic equation in the form x^2 + a*x + b = 0 and returns the
// real roots.
PolynomialRoots SolveQuadraticEquation(const double a, const double b) {
  // pre-computed constants
  double delta = a * a - 4.0 * b;
  // compute roots
  if (GetNearZero(delta)) {
    // two, identical real roots
    return {-0.5 * a, -0.5 * a};
  } else if (delta > 0) {
    // two real roots
    return {0.5 * (-a + sqrt(delta)), 0.5 * (-a - sqrt(delta))};
  }  // else only complex roots

  return {};
}

// Solves the Quadratic equation, but if a^2-4b < 0, will 0 that out and return
// real roots.
PolynomialRoots SafeSolveQuadraticEquation(const double a, const double b) {
  // pre-computed constants
  double delta = a * a - 4.0 * b;

  if (delta > 0) {
    // two real roots
    return {0.5 * (-a + sqrt(delta)), 0.5 * (-a - sqrt(delta))};
  }

  // two, identical real roots
  return {-0.5 * a, -0.5 * a};
}

// Solves the cubic equation in the form x^3 + a*x^2 + b*x + c = 0
// The return value is the number of real roots, and the roots themselves.
// If there are three real roots, all three values form the returned roots.
// If there is one real root and a pair of complex roots, then the real root is
// roots[0], and the pair of complex roots is given by roots[1] +/- i roots[2]
std::pair<int, PolynomialRoots> SolveCubicEquationWithComplexValues(
    const double a, const double b, const double c) {
  // pre-computed constants
  double delta_q = (Power2(a) - 3.0 * b) / 9.0;
  double delta_r = (a * (2.0 * Power2(a) - 9.0 * b) + 27.0 * c) / 54.0;

  // compute roots
  if (Power2(delta_r) < Power3(delta_q)) {
    // three real roots
    double phi = delta_r / sqrt(Power3(delta_q));
    if (phi < -1.0) {
      phi = -1.0;
    }
    if (phi > 1.0) {
      phi = 1.0;
    }
    phi = acos(phi);
    double a1 = a / 3.0;
    delta_q = -2.0 * sqrt(delta_q);
    return {3,
            {delta_q * cos(phi / 3.0) - a1,
             delta_q * cos((phi + 2.0 * M_PI) / 3.0) - a1,
             delta_q * cos((phi - 2.0 * M_PI) / 3.0) - a1}};
  }

  double delta_a =
      -std::cbrt(fabs(delta_r) + sqrt(Power2(delta_r) - Power3(delta_q)));
  if (delta_r < 0.0) {
    delta_a *= -1.0;
  }
  double delta_b = 0.0;
  if (delta_a != 0.0) {
    delta_b = delta_q / delta_a;
  }
  double delta_c = 0.5 * sqrt(3.0) * (delta_a - delta_b);
  double ceoff_a1 = a / 3.0;

  PolynomialRoots roots(3);
  roots[0] = delta_a + delta_b - ceoff_a1;
  roots[1] = -0.5 * (delta_a + delta_b) - ceoff_a1;
  if (fabs(delta_c) < kCubicDoubleRootsEpsilonScaleFactor *
                          std::numeric_limits<double>::epsilon()) {
    // two real roots
    roots[2] = roots[1];
    return {2, roots};
  }

  // one real roots x
  // two complex roots y +/- iz
  // roots[0] = x, roots[1] = y, roots[2] = z
  roots[2] = delta_c;
  return {1, roots};
}

// Solves the cubic equation in the form x^3 + a*x^2 + b*x + c = 0 and returns
// the real roots. Can return 1 or 3 real roots.
PolynomialRoots SolveCubicEquation(const double a, const double b,
                                   const double c) {
  auto [num_real_roots, roots] = SolveCubicEquationWithComplexValues(a, b, c);

  // Drop the complex stuff
  roots.resize(num_real_roots);
  return roots;
}

// Solves the quartic equation of the form x^4 + b*x^2 + d = 0
// Can return 4, 2, or zero real roots.
PolynomialRoots SolveBiQuadraticEquation(const double b, const double d) {
  // pre-computed constants
  double delta = b * b - 4.0 * d;
  if (fabs(delta) < kZeroThreshold && delta < 0.0) {
    delta = 0.0;
  }
  // compute roots
  if (delta >= 0) {
    // real roots
    double sqrt_delta = sqrt(delta);
    // two real roots for x^2
    double x2_root1, x2_root2;
    x2_root1 = 0.5 * (-b + sqrt_delta);
    x2_root2 = 0.5 * (-b - sqrt_delta);
    if (x2_root2 >= 0.0) {
      // four real roots
      double x_root1 = sqrt(x2_root1);
      double x_root2 = sqrt(x2_root2);
      return {-x_root1, x_root1, -x_root2, x_root2};
    } else if (x2_root1 >= 0.0) {
      // two real roots
      double x_root1 = sqrt(x2_root1);
      return {-x_root1, x_root1};
    }
  }
  // No real roots
  return {};
}

// Equation: x^4 + b*x^2 + c*x + d = 0
PolynomialRoots SolveDepressedQuarticEquation(const double b, const double c,
                                              const double d) {
  // check if the equation is biquadratic
  if (fabs(c) < kBiQuadraticScaleFactor * (fabs(b) + fabs(d))) {
    return SolveBiQuadraticEquation(b, d);
  }

  // solve resolvent cubic
  auto [resolvent_cubic_num_roots, resolvent_cubic_roots] =
      SolveCubicEquationWithComplexValues(2.0 * b, b * b - 4.0 * d, -c * c);
  // compute roots
  if (resolvent_cubic_num_roots > 1) {
    // sort the roots
    std::sort(resolvent_cubic_roots.begin(), resolvent_cubic_roots.end());
    if ((fabs(resolvent_cubic_roots[0]) < kZeroThreshold) &&
        (resolvent_cubic_roots[0] < 0.0)) {
      resolvent_cubic_roots[0] = 0.0;
    }
    if ((fabs(resolvent_cubic_roots[2]) < kZeroThreshold) &&
        (resolvent_cubic_roots[2] < 0.0)) {
      resolvent_cubic_roots[2] = 0.0;
    }
    if (resolvent_cubic_roots[0] >= 0.0) {
      // four positive roots
      if ((fabs(resolvent_cubic_roots[1]) < kZeroThreshold) &&
          (resolvent_cubic_roots[1] < 0.0)) {
        resolvent_cubic_roots[1] = 0.0;
      }
      double sqrt_res_root0 = sqrt(resolvent_cubic_roots[0]);
      double sqrt_res_root1 = sqrt(resolvent_cubic_roots[1]);
      double sqrt_res_root2 = sqrt(resolvent_cubic_roots[2]);
      if (c > 0.0) {
        return {0.5 * (-sqrt_res_root0 - sqrt_res_root1 - sqrt_res_root2),
                0.5 * (-sqrt_res_root0 + sqrt_res_root1 + sqrt_res_root2),
                0.5 * (sqrt_res_root0 - sqrt_res_root1 + sqrt_res_root2),
                0.5 * (sqrt_res_root0 + sqrt_res_root1 - sqrt_res_root2)};
      }
      return {0.5 * (-sqrt_res_root0 - sqrt_res_root1 + sqrt_res_root2),
              0.5 * (-sqrt_res_root0 + sqrt_res_root1 - sqrt_res_root2),
              0.5 * (sqrt_res_root0 - sqrt_res_root1 - sqrt_res_root2),
              0.5 * (sqrt_res_root0 + sqrt_res_root1 + sqrt_res_root2)};
    }

    // two pairs of complex roots, no real roots
    return {};
  }

  // two real roots and a pair of complex roots
  if (resolvent_cubic_roots[0] < 0.0) {
    resolvent_cubic_roots[0] *= -1.0;
  }

  double sqrt_res_root0 = sqrt(resolvent_cubic_roots[0]);
  double sqrt_real_part =
      std::sqrt(std::complex<double>(resolvent_cubic_roots[1],
                                     resolvent_cubic_roots[2]))
          .real();
  if (c > 0.0) {
    return {-0.5 * sqrt_res_root0 - sqrt_real_part,
            -0.5 * sqrt_res_root0 + sqrt_real_part};
  }

  return {0.5 * sqrt_res_root0 - sqrt_real_part,
          0.5 * sqrt_res_root0 + sqrt_real_part};
}

// One step of Newton Method for x^4 + a*x^3 + b*x^2 + c*x + d
double UpdateRootOneStepNewton(const double a, const double b, const double c,
                               const double d, const double root) {
  double first_derivative_func_value =
      ((4.0 * root + 3.0 * a) * root + 2.0 * b) * root + c;
  if (fabs(first_derivative_func_value) < kZeroThreshold) {
    // return without updating
    return root;
  }

  double func_value = (((root + a) * root + b) * root + c) * root + d;

  // update root
  return root - (func_value / first_derivative_func_value);
}

// Solves the quartic equation of the form x^4 + a*x^3 + b*x^2 + c*x + d = 0
// Can return 4, 2, or zero real roots.
PolynomialRoots SolveQuarticEquation(const double a, const double b,
                                     const double c, const double d) {
  // compute depressed quartic equation coefficients x^4 + b1*x^2 + c1*x + d1 =
  // 0
  double b1 = b - 0.375 * Power2(a);
  double c1 = c + 0.5 * a * (0.25 * Power2(a) - b);
  double d1 = d + 0.25 * a * (0.25 * b * a - (3.0 / 64.0) * Power3(a) - c);
  // compute roots
  auto roots = SolveDepressedQuarticEquation(b1, c1, d1);
  if (roots.empty()) {
    return roots;
  }

  for (double& root : roots) {
    root -= 0.25 * a;
  }

  // update real roots using one step of Newton method
  for (double& root : roots) {
    root = UpdateRootOneStepNewton(a, b, c, d, root);
  }
  return roots;
}

}  // namespace

PolynomialRoots CalculatePolynomialRoots(
    const std::array<double, 2>& coefficients) {
  if (!GetNearZero(coefficients[1])) {
    // Linear solve
    return {-coefficients[0] / coefficients[1]};
  }

  // Zeroth degree, empty
  return {};
}

PolynomialRoots CalculatePolynomialRoots(
    const std::array<double, 3>& coefficients) {
  if (!GetNearZero(coefficients[2])) {
    return SolveQuadraticEquation(coefficients[1] / coefficients[2],
                                  coefficients[0] / coefficients[2]);
  }

  return CalculatePolynomialRoots(Head<2>(coefficients));
}

PolynomialRoots SafeCalculatePolynomialRoots(
    const std::array<double, 3>& coefficients) {
  if (!GetNearZero(coefficients[2])) {
    return SafeSolveQuadraticEquation(coefficients[1] / coefficients[2],
                                      coefficients[0] / coefficients[2]);
  }

  // Make sure there are always 2 roots.
  PolynomialRoots roots = CalculatePolynomialRoots(Head<2>(coefficients));
  roots.resize(2, 0.);

  return roots;
}

PolynomialRoots CalculatePolynomialRoots(
    const std::array<double, 4>& coefficients) {
  if (!GetNearZero(coefficients[3])) {
    return SolveCubicEquation(coefficients[2] / coefficients[3],
                              coefficients[1] / coefficients[3],
                              coefficients[0] / coefficients[3]);
  }

  return CalculatePolynomialRoots(Head<3>(coefficients));
}

PolynomialRoots CalculatePolynomialRoots(
    const std::array<double, 5>& coefficients) {
  if (!GetNearZero(coefficients[4])) {
    return SolveQuarticEquation(
        coefficients[3] / coefficients[4], coefficients[2] / coefficients[4],
        coefficients[1] / coefficients[4], coefficients[0] / coefficients[4]);
  }

  return CalculatePolynomialRoots(Head<4>(coefficients));
}

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic
