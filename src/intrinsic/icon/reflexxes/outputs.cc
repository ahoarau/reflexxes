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

#include "intrinsic/icon/reflexxes/outputs.h"

#include <cstdint>
#include <string>
#include <utility>

#include "absl/log/check.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_format.h"
#include "absl/time/time.h"
#include "intrinsic/icon/reflexxes/constants_external.h"
#include "intrinsic/icon/reflexxes/inputs.h"
#include "intrinsic/icon/reflexxes/internal/input_output_utility_functions.h"
#include "intrinsic/icon/reflexxes/motion_polynomials.h"
#include "intrinsic/icon/reflexxes/profile.h"
#include "intrinsic/icon/reflexxes/status.h"

namespace intrinsic {
namespace reflexxes {

void Outputs::DOF::Scale(const double scale_factor) {
  CHECK(scale_factor > 0.0);

  new_position *= scale_factor;
  new_velocity *= scale_factor;
  new_acceleration *= scale_factor;

  // position extrema
  max_position *= scale_factor;
  min_position *= scale_factor;

  // motion polynomials
  motion_polynomials.Scale(scale_factor);
}

std::string Outputs::DOF::SubStep::GetProfileTraceString() const {
  std::string result;
  for (int i = 0; i < applied_profile_trace_size; ++i) {
    absl::StrAppend(&result, GetProfileString(applied_profile_trace[i]));
    if (i < applied_profile_trace_size - 1) {
      absl::StrAppend(&result, "->");
    }
  }

  return result;
}

static std::pair<bool, uint8_t> DecodeDecisionTrace(const uint8_t value) {
  return std::pair(value & 128, value & 127);
}

std::string Outputs::DOF::SubStep::DecisionTrace::GetTraceString(
    int count) const {
  std::string result;
  for (int i = 0; i < count; ++i) {
    auto [decision, num] = DecodeDecisionTrace(decisions_[i]);
    absl::StrAppendFormat(&result, "%d%s", num, decision ? "y" : "n");
    if (i < count - 1) {
      absl::StrAppend(&result, "->");
    }
  }

  return result;
}

Status Outputs::ComputeNextStateOfMotion() {
  time_counter_++;

  // Use absl::Duration to avoid unnecessary rounding errors when computing the
  // new time.
  absl::Duration new_time = time_counter_ * absl::Seconds(GetCycleTime());

  if (new_time > absl::Seconds(kMaxMinExecutionTime)) {
    SetStatus(Status::kErrorUserTimeOutOfRange);
    return Status::kErrorUserTimeOutOfRange;
  }

  // compute new state
  for (DOF& dof : GetDOFs()) {
    MotionPolynomials::MotionState motion_state =
        dof.motion_polynomials.GetStateOfMotionAtTime(
            absl::ToDoubleSeconds(new_time));

    dof.new_position = motion_state.position;
    dof.new_velocity = motion_state.velocity;
    dof.new_acceleration = motion_state.acceleration;
  }

  if (new_time < absl::Seconds(GetSyncTime())) {
    SetStatus(Status::kWorking);
  } else {
    SetStatus(Status::kFinalStateReached);
  }

  // indicate that a valid output is available
  is_valid_output_available_ = true;
  return status_;
}

bool Outputs::CopyNewStateToCurrentState(Inputs& inputs) const {
  if (inputs.GetNumberOfDOFs() != GetNumberOfDOFs()) {
    return false;
  }

  internal::CopyNewStateToInputs(*this, inputs);
  return true;
}

Outputs::Outputs(const int num_dofs, const double cycle_time)
    : cycle_time_(cycle_time), dofs_(num_dofs) {
  for (int i = 0; i < num_dofs; i++) {
    dofs_[i].index = i;
  }
}

Outputs::Outputs(const MaxDOFFixedVector<DOF>& dofs, const double cycle_time,
                 const double sync_time, const uint64_t time_counter,
                 const Status status, const bool is_valid_output_available,
                 const bool phase_sync_enabled, const int phase_sync_dof_index)
    : cycle_time_(cycle_time),
      sync_time_(sync_time),
      time_counter_(time_counter),
      status_(status),
      is_valid_output_available_(is_valid_output_available),
      phase_sync_enabled_(phase_sync_enabled),
      phase_sync_dof_index_(phase_sync_dof_index),
      dofs_(dofs) {
  for (int i = 0; i < dofs.size(); i++) {
    dofs_[i].index = i;
  }
}

}  // namespace reflexxes
}  // namespace intrinsic
