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

#ifndef INTRINSIC_ICON_REFLEXXES_STATUS_H_
#define INTRINSIC_ICON_REFLEXXES_STATUS_H_

namespace intrinsic {
namespace reflexxes {

// Status codes for the result of the algorithm.
enum class Status {
  kWorking = 0,
  kFinalStateReached = 1,
  kErrorNumberOfDofs = -101,
  kErrorCycleTime = -102,
  kErrorInvalidInputValues = -103,
  kErrorInvalidTargetState = -104,
  kErrorInvalidConstraints = -105,
  kErrorInvalidScaleOfInputValues = -106,
  kErrorExecutionTimeCalculation = -107,
  kErrorSynchronization = -108,
  kErrorNoPhaseSynchronization = -109,
  kErrorExecutionTimeTooBig = -110,
  kErrorUserTimeOutOfRange = -111,
  kErrorPositionalLimits = -112,
  kErrorUndefined = -999,
};

// Get a C-string containing the status message.
constexpr const char* GetStatusString(const Status status) {
  switch (status) {
    case Status::kWorking:
      return "kWorking: The Online Trajectory Generation algorithm is working; "
             "the final state of motion has not been reached yet.";
    case Status::kFinalStateReached:
      return "kFinalStateReached: The desired final state of motion has been "
             "reached.";
    case Status::kErrorNumberOfDofs:
      return "kErrorNumberOfDofs: The number of degree of freedom of the "
             "input parameters, the output parameters, and the Online "
             "Trajectory Generation algorithm do not match.";
    case Status::kErrorCycleTime:
      return "kErrorCycleTime: The cycle time of the input parameters, the "
             "output parameters, and the Online Trajectory Generation "
             "algorithm do not match.";
    case Status::kErrorInvalidInputValues:
      return "kErrorInvalidInputValues: The applied input values are "
             "invalid.";
    case Status::kErrorInvalidTargetState:
      return "kErrorInvalidTargetState: Target State in Input Parameters is "
             "invalid";
    case Status::kErrorInvalidConstraints:
      return "kErrorInvalidConstraints: Constraints in Input Parameters are "
             "invalid";
    case Status::kErrorInvalidScaleOfInputValues:
      return "kErrorInvalidScaleOfInputValues: Scale of Input Values is "
             "too big";
    case Status::kErrorExecutionTimeCalculation:
      return "kErrorExecutionTimeCalculation (STEP 1): An error occurred "
             "during the calculation of the synchronization time.";
    case Status::kErrorSynchronization:
      return "kErrorSynchronization (STEP 2): An error occurred during the "
             "synchronization of the trajectory.";
    case Status::kErrorNoPhaseSynchronization:
      return "kErrorNoPhaseSynchronization: The input flag "
             "Flags::kOnlyPhaseSynchronization is set, but a "
             "phase-synchronized trajectory cannot be executed.";
    case Status::kErrorExecutionTimeTooBig:
      return "kErrorExecutionTimeTooBig: The execution time of the computed "
             "trajectory is too big (>kMaxAllowedSyncTimeSeconds seconds).";
    case Status::kErrorUserTimeOutOfRange:
      return "kErrorUserTimeOutOfRange: The sample time for the previously "
             "computed trajectory is out of range.";
    case Status::kErrorPositionalLimits:
      return "kErrorPositionalLimits: The computed trajectory will exceed the "
             "positional limits.";
    case Status::kErrorUndefined:
      return "kErrorUndefined: An unknown error has occurred.";
  }
}

}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_STATUS_H_
