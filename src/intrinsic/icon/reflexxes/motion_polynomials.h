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

#ifndef INTRINSIC_ICON_REFLEXXES_MOTION_POLYNOMIALS_H_
#define INTRINSIC_ICON_REFLEXXES_MOTION_POLYNOMIALS_H_

#include <cmath>
#include <cstdio>
#include <limits>

#include "absl/types/span.h"
#include "intrinsic/icon/reflexxes/constants_external.h"
#include "intrinsic/util/fixed_vector.h"

namespace intrinsic {
namespace reflexxes {

// Represents a motion trajectory comprised of a set of time continuous
// segments, where each segment defines a different motion profile.
class MotionPolynomials {
 public:
  struct Segment {
    double position = 0.;
    double velocity = 0.;
    double acceleration = 0.;
    double jerk = 0.;

    // Starting time for the segment in seconds.
    double start_time = 0.0;

    // Ending time for the segment in seconds.
    double end_time = 0.0;
  };

  // Struct describing motion state of a single DOF containing its position,
  // velocity and acceleration values.
  struct MotionState {
    double position = 0.0;
    double velocity = 0.0;
    double acceleration = 0.0;
  };

  // Struct describing the extrema of a single category of a single DOF and the
  // times at which they occur.
  struct Extrema {
    double max_value = -std::numeric_limits<double>::max();
    double min_value = std::numeric_limits<double>::max();
    double max_time = 0.0;
    double min_time = 0.0;
  };

  absl::Span<Segment> GetSegments() { return absl::MakeSpan(segments_); }

  absl::Span<const Segment> GetSegments() const { return segments_; }

  int GetSegmentCount() const { return segments_.size(); }

  bool HasSegments() const { return GetSegmentCount() > 0; }

  // Sets the number of segments to count.  If count is less than the existing
  // number of segments, N, then N - count segments are removed from the end of
  // the list.  If the count is equal to or greater than the number of segments,
  // no action is taken.
  void TrimSegmentCount(int count);

  // Adds a segment to the end of the trajectory.
  void AddSegment(const Segment& segment);

  // Scales the internal segment PVAJ values by a constant amount.
  void Scale(double scale_factor);

  // Time scales the segments by adjusting each segment's start and end time by
  // the given time_scale factor, and scaling the PVAJ attributes of the
  // segments to match accordingly.
  void TimeScale(double time_scale);

  // Scales segments such that the first segment is set to the given PVA values,
  // and the remaining segments are scaled by the given scale factor.
  void ScaleWithNewInitialValues(double scale_factor, double initial_position,
                                 double initial_velocity,
                                 double initial_acceleration);

  // Returns the state of motion (position, velocity and acceleration) of a
  // single degree of freedom at a given time.  If the given time is not
  // encapsulated in the range of times described by this trajectory, a default
  // constructed MotionState is returned instead.
  MotionState GetStateOfMotionAtTime(double time_in_seconds) const;

  // Returns an Extrema object that contains the maximum and minimum
  // position values and their associated times among all segments of the
  // position polynomial
  Extrema GetPositionExtrema() const {
    return GetPositionExtrema(GetSegmentCount());
  }

  // Returns an Extrema object that contains the maximum and minimum
  // position values and their associated times up to num_segments of the
  // position polynomial
  Extrema GetPositionExtrema(int num_segments) const;

  // Returns an Extrema object that contains the maximum and minimum
  // velocity values and their associated times for all segments of the velocity
  // polynomial
  Extrema GetVelocityExtrema() const;

  // Returns an Extrema object that contains the maximum and minimum
  // acceleration values and their associated times for all segments of the
  // acceleration polynomial
  Extrema GetAccelerationExtrema() const;

  // Returns an Extrema object that contains the maximum and minimum jerk
  // values from the acceleration polynomial
  Extrema GetJerkExtrema() const;

 private:
  FixedVector<Segment, kMaxNumPolynomials> segments_;
};

}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_MOTION_POLYNOMIALS_H_
