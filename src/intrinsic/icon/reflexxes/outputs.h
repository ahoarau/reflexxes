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

#ifndef INTRINSIC_ICON_REFLEXXES_OUTPUTS_H_
#define INTRINSIC_ICON_REFLEXXES_OUTPUTS_H_

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "absl/log/check.h"
#include "absl/types/span.h"
#include "intrinsic/icon/reflexxes/constants_external.h"
#include "intrinsic/icon/reflexxes/inputs.h"
#include "intrinsic/icon/reflexxes/motion_polynomials.h"
#include "intrinsic/icon/reflexxes/profile.h"
#include "intrinsic/icon/reflexxes/status.h"

namespace intrinsic {
namespace reflexxes {

// Base class for the results of the Online Trajectory Generation algorithm.
// This class cannot be constructed directly, instead use one of the child
// classes, PositionOutputs or VelocityOutputs.
class Outputs {
 public:
  // Structure containing the output parameters for a single DOF
  struct DOF {
    // Structure containing the substep output parameters for a single DOF
    // This represents the result of a single substep of execution (e.g.
    // Position Step 1B, Velocity Step 2).
    struct SubStep {
      // Helper class for tracing decisions in a decision tree.
      class DecisionTrace {
       public:
        // Encodes a decision trace at the specified index.
        void AddDecisionTrace(int index, const bool decision,
                              const uint8_t num) {
          CHECK(index < sizeof(decisions_));
          decisions_[index] =
              static_cast<int8_t>((decision ? 128 : 0) | (num & 127));
        }

        // Gets the string representing the decision trace stored in this
        // substep.  Gets 'size' number of decisions.
        std::string GetTraceString(int count) const;

       private:
        // The character array of the sequence of decision nodes visited.
        uint8_t decisions_[kTraceSize];
      };

      // The profile id applied.
      Profile applied_profile = Profile::kUndefined;

      // The result (success or failure).
      bool result = false;

      // The character array of the sequence of applied profiles.
      using ProfileTrace = std::array<Profile, kTraceSize>;
      ProfileTrace applied_profile_trace;

      // Used size of the applied_provile_trace array.
      int8_t applied_profile_trace_size = 0;

      // The trace of decisions for this substep.
      DecisionTrace decision_trace;
      // Used size of the decision_trace_size array.
      int8_t decision_trace_size = 0;

      // Execution time for this step.
      double step_exec_time = 0.;

      // Gets the string representing the profile trace stored in this substep.
      std::string GetProfileTraceString() const;
    };

    // The index of this object in the larger structure.
    int index = 0;

    // Whether or not this dof was selected for the output.
    bool selected = false;

    double new_position = 0.0;
    double new_velocity = 0.0;
    double new_acceleration = 0.0;

    double min_execution_time = 0.0;

    double inoperative_begin_execution_time = 0.0;

    double inoperative_end_execution_time = 0.0;

    // Execution time of the degree of freedom in the case non-synchronized
    // trajectories.
    double execution_time = 0.0;

    // The motion polynomials object which contains the vector of polynomial
    // coefficients of position, velocity and acceleration.
    MotionPolynomials motion_polynomials;

    // Scale factor used for phase synchronization.
    double phase_sync_scale = 1.0;

    // Time at which the degree of freedom reaches its minimum position during
    // the execution of the calculated trajectory.
    double min_position_extrema_time = 0.0;

    // Time at which the degree of freedom reaches its maximum position during
    // the execution of the calculated trajectory.
    double max_position_extrema_time = 0.0;

    // Minimum position for the degree of freedom that occurs during the
    // execution of the calculated trajectory.
    double min_position = 0.0;

    // Maximum position for the degree of freedom that occurs during the
    // execution of the calculated trajectory.
    double max_position = 0.0;

    // Time at which the degree of freedom reaches its minimum velocity during
    // the execution of the calculated trajectory.
    double min_velocity_extrema_time = 0.0;

    // Time at which the degree of freedom reaches its maximum velocity during
    // the execution of the calculated trajectory.
    double max_velocity_extrema_time = 0.0;

    // Minimum velocity for the degree of freedom that occurs during the
    // execution of the calculated trajectory.
    double min_velocity = 0.0;

    // Maximum velocity for the degree of freedom that occurs during the
    // execution of the calculated trajectory.
    double max_velocity = 0.0;

    // Time at which the degree of freedom reaches its minimum acceleration
    // during the execution of the calculated trajectory.
    double min_acceleration_extrema_time = 0.0;

    // Time at which the degree of freedom reaches its maximum acceleration
    // during the execution of the calculated trajectory.
    double max_acceleration_extrema_time = 0.0;

    // Minimum acceleration for the degree of freedom that occurs during the
    // execution of the calculated trajectory.
    double min_acceleration = 0.0;

    // Maximum acceleration for the degree of freedom that occurs during the
    // execution of the calculated trajectory.
    double max_acceleration = 0.0;

    // Time at which the degree of freedom reaches its minimum jerk
    // during the execution of the calculated trajectory.
    double min_jerk_extrema_time = 0.0;

    // Time at which the degree of freedom reaches its maximum jerk
    // during the execution of the calculated trajectory.
    double max_jerk_extrema_time = 0.0;

    // Minimum jerk for the degree of freedom that occurs during the
    // execution of the calculated trajectory.
    double min_jerk = 0.0;

    // Maximum jerk for the degree of freedom that occurs during the
    // execution of the calculated trajectory.
    double max_jerk = 0.0;

    // Substep output data for Step1/Step1A (velocity vs. position).
    SubStep step1a;

    // Substep output data for Step1B (velocity vs. position).
    SubStep step1b;

    // Substep output data for Step1C (velocity vs. position).
    SubStep step1c;

    // Substep output data for Step2.
    SubStep step2;

    // Returns true if the exec_time is within the inoperative window of this
    // dof.
    bool IsWithinAnInoperativeTimeInterval(const double exec_time) const {
      return (exec_time > inoperative_begin_execution_time &&
              exec_time < inoperative_end_execution_time);
    }

    // Sets the execution time if this dof is selected.
    void SetExecutionTimeIfSelected(const double value) {
      execution_time = selected ? value : 0.0;
    }

    void Scale(double scale_factor);
  };

  // Returns the array of all dofs for this output.
  absl::Span<DOF> GetDOFs() { return absl::MakeSpan(dofs_); }

  // Returns the array of all dofs for this output.
  absl::Span<const DOF> GetDOFs() const { return absl::MakeConstSpan(dofs_); }

  // Computes the next state of motion for all degrees of freedom.
  // This method implements Step3 of the Reflexxes Motion Libraries 2.0.
  Status ComputeNextStateOfMotion();

  // Copy the new state from output parameters to the current state in the
  // given input parameters.
  // Returns false if the number of DOFs do not match.
  bool CopyNewStateToCurrentState(Inputs& inputs) const;

  // Returns true if exec_time is within the inoperative time window of any
  // selected dof.
  bool IsWithinAnInoperativeTimeInterval(const double exec_time) const {
    for (const auto& dof_output : dofs_) {
      if (dof_output.selected &&
          dof_output.IsWithinAnInoperativeTimeInterval(exec_time))
        return true;
    }

    return false;
  }

  // Sets the execution time of any selected dof to value, and non selected
  // dofs to 0.
  void SetExecutionTimeIfSelected(const double value) {
    for (auto& dof_output : GetDOFs()) {
      dof_output.SetExecutionTimeIfSelected(value);
    }
  }

  // Get the number of degrees of freedom.
  int GetNumberOfDOFs() const { return dofs_.size(); }

  // Get the cycle time.
  double GetCycleTime() const { return cycle_time_; }

  // Get the synchronization time.
  double GetSyncTime() const { return sync_time_; }

  // Set the synchronization time.
  void SetSyncTime(const double sync_time) { sync_time_ = sync_time; }

  // Get the status of execution.
  Status GetStatus() const { return status_; }

  // Set the status of execution.
  void SetStatus(const Status& status) { status_ = status; }

  // Returns true if this output was previously set via a valid compute call.
  bool IsValidOutputAvailable() const { return is_valid_output_available_; }

  // Get the phase synchronization boolean flag.
  bool IsPhaseSyncEnabled() const { return phase_sync_enabled_; }

  // Sets whether or not the phase sync is on, and if so, which dof to use.
  void EnablePhaseSync(const int index) {
    phase_sync_enabled_ = true;
    phase_sync_dof_index_ = index;
  }

  // Set the phase synchronization boolean flag.
  void DisablePhaseSync() { phase_sync_enabled_ = false; }

  // Get the DOF used for phase synchronization if its on.
  int GetPhaseSyncDOFIndex() const { return phase_sync_dof_index_; }

  // Counter for ticks in cycle time, updated with each call to
  // computeNextStateOfMotion.
  uint64_t GetTimeCounter() const { return time_counter_; }

 protected:
  // Constructor of class Outputs
  Outputs(int num_dofs, double cycle_time);

  // Constructor where everything is set.
  Outputs(const MaxDOFFixedVector<DOF>& dofs, double cycle_time,
          double sync_time, uint64_t time_counter, Status status,
          bool is_valid_output_available, bool phase_sync_enabled,
          int phase_sync_dof_index);

 private:
  double cycle_time_;

  // The synchronization time in seconds. It contains the maximum execution
  // time when the trajectories are asynchronous
  double sync_time_ = 0.0;

  // Counter for the time at which the new state is calculated
  uint64_t time_counter_ = 0;

  // Reports status of the Reflexxes Motion Library
  Status status_ = Status::kErrorUndefined;

  // Boolean flag that indicates whether a valid output is available
  bool is_valid_output_available_ = false;

  // True if the motion trajectories are phase synchronized
  bool phase_sync_enabled_ = false;

  // DOF used for phase synchronization
  int phase_sync_dof_index_ = 0;

  // Vector of DOF objects that contains the output parameters for all degrees
  // of freedom
  MaxDOFFixedVector<DOF> dofs_;
};

}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_OUTPUTS_H_
