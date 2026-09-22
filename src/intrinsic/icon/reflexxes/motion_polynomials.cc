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

#include "intrinsic/icon/reflexxes/motion_polynomials.h"

#include <array>
#include <functional>

#include "absl/algorithm/container.h"
#include "absl/log/check.h"
#include "absl/types/span.h"
#include "intrinsic/icon/reflexxes/constants_external.h"
#include "intrinsic/icon/reflexxes/internal/functional.h"
#include "intrinsic/icon/reflexxes/internal/math.h"
#include "intrinsic/icon/reflexxes/internal/polynomial_solvers.h"

namespace intrinsic {
namespace reflexxes {
namespace {
// Bring in some common types and functions.
using Extrema = MotionPolynomials::Extrema;
using Segment = MotionPolynomials::Segment;
using internal::EvalPolynomial;
using internal::Power2;
using internal::Power3;

std::array<double, 4> CalculatePositionCoeffs(const Segment& segment) {
  return {segment.position, segment.velocity, segment.acceleration / 2.0,
          segment.jerk / 6.0};
}

std::array<double, 4> CalculateVelocityCoeffs(const Segment& segment) {
  return {segment.velocity, segment.acceleration, segment.jerk / 2.0, 0};
}

std::array<double, 4> CalculateAccelerationCoeffs(const Segment& segment) {
  return {segment.acceleration, segment.jerk, 0, 0};
}

double CalculatePosition(const Segment& segment, const double t) {
  return EvalPolynomial(CalculatePositionCoeffs(segment),
                        t - segment.start_time);
}

double CalculateVelocity(const Segment& segment, const double t) {
  return EvalPolynomial(CalculateVelocityCoeffs(segment),
                        t - segment.start_time);
}

double CalculateAcceleration(const Segment& segment, const double t) {
  return EvalPolynomial(CalculateAccelerationCoeffs(segment),
                        t - segment.start_time);
}
}  // namespace

void MotionPolynomials::Scale(double scale_factor) {
  for (Segment& segment : segments_) {
    segment.position *= scale_factor;
    segment.velocity *= scale_factor;
    segment.acceleration *= scale_factor;
    segment.jerk *= scale_factor;
  }
}

void MotionPolynomials::TimeScale(double time_scale) {
  CHECK(time_scale > 0.0);

  if (!HasSegments()) {
    return;
  }

  // Scale the first one
  Segment& first_segment = segments_.front();
  first_segment.start_time *= time_scale;
  first_segment.end_time *= time_scale;
  first_segment.velocity /= time_scale;
  first_segment.acceleration /= Power2(time_scale);
  first_segment.jerk /= Power3(time_scale);

  // Iterate through the rest of the segments and their previous segments and
  // calculate their PVA based off the previous segment.
  for (auto [previous_segment, segment] : internal::Zip(
           GetSegments(), GetSegments().last(GetSegmentCount() - 1))) {
    segment.start_time *= time_scale;
    segment.end_time *= time_scale;
    segment.position = CalculatePosition(previous_segment, segment.start_time);
    segment.velocity = CalculateVelocity(previous_segment, segment.start_time);
    segment.acceleration =
        CalculateAcceleration(previous_segment, segment.start_time);
    segment.jerk /= Power3(time_scale);
  }
}

void MotionPolynomials::ScaleWithNewInitialValues(double scale_factor,
                                                  double initial_position,
                                                  double initial_velocity,
                                                  double initial_acceleration) {
  if (!HasSegments()) {
    return;
  }

  // Set the initial values of the trajectory and scale the first jerk.
  Segment& first_segment = segments_.front();
  first_segment.position = initial_position;
  first_segment.velocity = initial_velocity;
  first_segment.acceleration = initial_acceleration;
  first_segment.jerk *= scale_factor;

  // Iterate through the rest of the segments and their previous segments and
  // calculate their PVA based off the previous segment.
  for (auto [previous_segment, segment] : internal::Zip(
           GetSegments(), GetSegments().last(GetSegmentCount() - 1))) {
    segment.position = CalculatePosition(previous_segment, segment.start_time);
    segment.velocity = CalculateVelocity(previous_segment, segment.start_time);
    segment.acceleration =
        CalculateAcceleration(previous_segment, segment.start_time);
    segment.jerk *= scale_factor;
  }
}

void MotionPolynomials::TrimSegmentCount(int count) {
  if (count >= 0 && count < segments_.size()) {
    segments_.resize(count);
  }
}

void MotionPolynomials::AddSegment(const Segment& segment) {
  CHECK_LT(GetSegmentCount(), kMaxNumPolynomials);
  segments_.push_back(segment);
}

MotionPolynomials::MotionState MotionPolynomials::GetStateOfMotionAtTime(
    double time_in_seconds) const {
  MotionState motion_state;

  // find the segment for the time requested
  const auto iter =
      absl::c_find_if(segments_, [time_in_seconds](const Segment& segment) {
        return time_in_seconds >= segment.start_time &&
               time_in_seconds < segment.end_time;
      });

  // compute the state
  if (iter != segments_.end()) {
    motion_state.position = CalculatePosition(*iter, time_in_seconds);
    motion_state.velocity = CalculateVelocity(*iter, time_in_seconds);
    motion_state.acceleration = CalculateAcceleration(*iter, time_in_seconds);
  }

  return motion_state;
}

namespace {
using internal::Head;
using internal::PolynomialRoots;

// Calculates the roots, and then offsets the roots result by the start time of
// the segment.
PolynomialRoots CalculatePolynomialRoots(
    const std::array<double, 3>& coefficients, const double start_time) {
  PolynomialRoots roots = internal::CalculatePolynomialRoots(coefficients);
  for (double& root : roots) {
    root += start_time;
  }
  return roots;
}

PolynomialRoots CalculateVelocityRoots(const Segment& segment) {
  return CalculatePolynomialRoots(Head<3>(CalculateVelocityCoeffs(segment)),
                                  segment.start_time);
}

PolynomialRoots CalculateAccelerationRoots(const Segment& segment) {
  return CalculatePolynomialRoots(Head<3>(CalculateAccelerationCoeffs(segment)),
                                  segment.start_time);
}

void MaybeUpdateExtrema(const double value, const double time,
                        Extrema& to_update) {
  if (value > to_update.max_value) {
    to_update.max_value = value;
    to_update.max_time = time;
  }
  if (value < to_update.min_value) {
    to_update.min_value = value;
    to_update.min_time = time;
  }
}

// Common function for getting the extrema of a property of the trajectory.
// Calculates the extrema over the set of passed in segments using the
// calculate_value function to get the value of the property for any time T from
// any segment.
// calculate_roots is an optional function pointer which may or may not be
// passed in. (e.g. getting the acceleration extrema doesn't use roots)
Extrema GetExtrema(
    absl::Span<const Segment> segments,
    const std::function<double(const Segment&, double)>& calculate_value,
    const std::function<PolynomialRoots(const Segment&)>& calculate_roots) {
  if (segments.empty()) {
    return Extrema{};
  }

  Extrema result;

  // initialize the extrema to the current position, i.e. position at zero
  // time
  double start_time = segments[0].start_time;
  double start_value = result.min_value = result.max_value =
      calculate_value(segments[0], start_time);

  // iterate through all polynomial segments and find the min and max values
  for (const Segment& segment : segments) {
    // compute start and end positions for each polynomial segment
    double end_time = segment.end_time;
    double end_value = calculate_value(segment, end_time);

    // find the max and min among start and end values
    MaybeUpdateExtrema(start_value, start_time, result);
    MaybeUpdateExtrema(end_value, end_time, result);

    if (calculate_roots != nullptr) {
      // compute extrema time as roots of velocity polynomial
      PolynomialRoots roots = calculate_roots(segment);

      // check if valid roots are found
      for (double root : roots) {
        if ((root >= start_time) && (root <= end_time)) {
          MaybeUpdateExtrema(calculate_value(segment, root), root, result);
        }
      }
    }

    // Reset start value/time for the next loop
    start_time = end_time;
    start_value = end_value;
  }

  return result;
}

}  // namespace

Extrema MotionPolynomials::GetPositionExtrema(const int num_segments) const {
  return GetExtrema(GetSegments().first(num_segments), &CalculatePosition,
                    &CalculateVelocityRoots);
}

Extrema MotionPolynomials::GetVelocityExtrema() const {
  return GetExtrema(GetSegments(), &CalculateVelocity,
                    &CalculateAccelerationRoots);
}

Extrema MotionPolynomials::GetAccelerationExtrema() const {
  return GetExtrema(GetSegments(), &CalculateAcceleration, nullptr);
}

Extrema MotionPolynomials::GetJerkExtrema() const {
  Extrema result;
  // iterate through all polynomial segments and find the min and max values
  for (const Segment& segment : segments_) {
    MaybeUpdateExtrema(segment.jerk, segment.start_time, result);
  }

  return result;
}

}  // namespace reflexxes
}  // namespace intrinsic
