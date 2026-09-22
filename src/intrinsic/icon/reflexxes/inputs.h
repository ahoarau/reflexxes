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

#ifndef INTRINSIC_ICON_REFLEXXES_INPUTS_H_
#define INTRINSIC_ICON_REFLEXXES_INPUTS_H_

#include <math.h>

#include <utility>
#include <vector>

#include "absl/types/span.h"
#include "intrinsic/icon/reflexxes/constants_external.h"
#include "intrinsic/util/fixed_vector.h"

namespace intrinsic {
namespace reflexxes {

// An array class of type T with enough slots for the max number of dofs.
template <typename T>
using MaxDOFFixedVector = FixedVector<T, kMaxDofs>;

// Base class for the inputs of the Online Trajectory Generation algorithm.
// This class cannot be constructed directly, instead use one of the child
// classes, PositionInputs or VelocityInputs.
class Inputs {
 public:
  virtual ~Inputs() = default;

  // Error codes for invalid input values.
  enum class ErrorCodeForInvalidInputValues {
    // No error has occurred.
    kNoError = 0,
    // The value for the maximum velocity is not greater than zero.
    kMaxVelocity = 1,
    // The value for the minimum velocity is not less than zero.
    kMinVelocity = 2,
    // The value for the maximum acceleration is not greater than zero.
    kMaxAcceleration = 3,
    // The value for the minimum acceleration is not less than zero.
    kMinAcceleration = 4,
    // The value for the maximum jerk is not greater than zero.
    kMaxJerk = 5,
    // The value for the minimum jerk is not less than zero.
    kMinJerk = 6,
    // The value for the target velocity is outside the range of the min/max
    // velocity.
    kTargetVelocity = 7,
    // The ratio of the absolute value of the largest constraint limit
    // (velocity, acceleration, jerk) to the smallest constraint limit is
    // greater than kPositionValidityMagnitude or kPositionValidityMagnitude for
    // position and velocity computations, respectively.  This will create
    // numerical instabilities and a result cannot be guaranteed.
    kOrderOfMagnitude = 8,
    // The value for the user-specified synchronization time is too large.
    kMinimumSyncTime = 9,
    // The value for the positional limits are too large in magnitude.
    kPositionLimitTooLarge = 10,
    // The value for the jerk limits are too large in magnitude.
    kJerkLimitTooLarge = 11,
  };

  // Structure containing the input parameters for a single DOF.
  struct DOF {
    // The index of this object in the larger structure.
    int index = 0;

    // Whether or not this DOF should be considered during execution of the
    // algorithm.  Generally all DOFs passed into the Compute calls are
    // selected, but in some cases this may not be desired, or in others the
    // algorithm (under specific runtime flags) may automatically deselect the
    // DOF.
    bool selected = true;

    double position = 0.0;
    double velocity = 0.0;
    double acceleration = 0.0;

    double target_position = 0.0;
    double target_velocity = 0.0;

    // The target velocity is used when fallback strategy is executed.
    double alt_target_velocity = 0.0;

    double max_position = 0.0;
    double min_position = 0.0;
    double max_velocity = 0.0;
    double min_velocity = 0.0;
    double max_acceleration = 0.0;
    double min_acceleration = 0.0;
    double max_jerk = 0.0;
    double min_jerk = 0.0;

    // Scales the DOF to the given scale factor.
    void ScaleValues(double scale_factor);
  };

  // Get all the DOFs.
  absl::Span<const DOF> GetDOFs() const { return absl::MakeSpan(dofs_); }

  // Get all the DOFs.
  absl::Span<DOF> GetDOFs() { return absl::MakeSpan(dofs_); }

  // Checks validity of the input.
  //
  // This is an important method to ensure numerical robustness of the
  // Online Trajectory Generation algorithm. Only if the result of this
  // method is true, a correct computation of output values can be
  // obtained. If the result is false, the Online Trajectory Generation
  // Algorithm will try compute correct output values, and in many cases,
  // this is possible, but it is not guaranteed. All input values have to
  // be within a proper order of magnitude (cf. kPositionValidityMagnitude),
  // the kinematic motion constraints, that is, the maximum values for
  // velocity, acceleration, and jerk, have to be positive, the minimum values
  // for velocity, acceleration, and jerk, have to be negative, and the target
  // velocities have to be within their corresponding minimum and maximum
  // bounds.
  //
  // Returns any error code or Inputs::ErrorCodeForInvalidInputValues::kNoError,
  // and the index of the DOF causing the error.  If there is no error or the
  // error is not specific to a DOF the second argument of the return value will
  // be -1.
  virtual std::pair<ErrorCodeForInvalidInputValues, int> CheckForValidity()
      const = 0;

  // Returns true if CheckForValidity() returns
  // Inputs::ErrorCodeForInvalidInputValues::kNoError.
  bool IsValid() const {
    return CheckForValidity().first == ErrorCodeForInvalidInputValues::kNoError;
  }

  // Returns true if the scale of the inputs is incorrect.  (CheckForValidity
  // returns kOrderOfMagnitude or kMinimumSyncTime).
  bool IsScaleOfInputsValid() const {
    ErrorCodeForInvalidInputValues code = CheckForValidity().first;
    return code != ErrorCodeForInvalidInputValues::kOrderOfMagnitude &&
           code != ErrorCodeForInvalidInputValues::kMinimumSyncTime;
  }

  // Sets the optional parameter for the minimum synchronization time in seconds
  void SetMinimumSynchronizationTime(const double sync_time) {
    min_sync_time_ = sync_time;
  }

  // Get the optional parameter for the minimum synchronization time in seconds.
  double GetMinimumSynchronizationTime() const { return min_sync_time_; }

  // Get the number of degrees of freedom.
  int GetNumberOfDOFs() const { return dofs_.size(); }

  // Get the cycle time.
  double GetCycleTime() const { return cycle_time_; }

 protected:
  // Constructor of class Inputs.
  Inputs(int num_dofs, double cycle_time);

  // Constructor of class Inputs
  Inputs(const MaxDOFFixedVector<DOF>& dofs, double cycle_time,
         double min_sync_time);

 private:
  // Vector of DOF object that contains the input parameters for all degrees of
  // freedom
  MaxDOFFixedVector<DOF> dofs_;

  double cycle_time_;

  // Optional minimum execution time in seconds specified by the user
  double min_sync_time_;
};

// flip signs of current position, velocity and acceleration, flip and swap
// limits for velocity, acceleration, and jerk
Inputs::DOF FlipInputParameters(const Inputs::DOF& i);

}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_INPUTS_H_
