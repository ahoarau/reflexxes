// Minimal 3-DOF position-based Reflexxes Type IV example.

#include <cstdio>

#include "intrinsic/icon/reflexxes/reflexxes_api.h"

namespace rx = intrinsic::reflexxes;

int main() {
  constexpr int kNumDofs = 3;
  constexpr double kCycleTime = 0.001;

  rx::State state(kNumDofs, kCycleTime);
  rx::PositionInputs inputs(kNumDofs, kCycleTime);
  rx::PositionOutputs outputs(kNumDofs, kCycleTime);
  rx::PositionFlags flags;

  const double targets[kNumDofs] = {1.0, -0.5, 2.0};
  for (int i = 0; i < kNumDofs; ++i) {
    rx::Inputs::DOF& dof = inputs.GetDOFs()[i];
    dof.index = i;
    dof.target_position = targets[i];
    dof.max_position = 10.0;
    dof.min_position = -10.0;
    dof.max_velocity = 1.0;
    dof.min_velocity = -1.0;
    dof.max_acceleration = 2.0;
    dof.min_acceleration = -2.0;
    dof.max_jerk = 10.0;
    dof.min_jerk = -10.0;
  }

  rx::Status status = rx::Status::kWorking;
  int cycles = 0;
  while (status == rx::Status::kWorking) {
    status = rx::ComputePosition(inputs, flags, outputs, state);
    if (status != rx::Status::kWorking &&
        status != rx::Status::kFinalStateReached) {
      std::fprintf(stderr, "Error: %s\n", rx::GetStatusString(status));
      return 1;
    }
    if (cycles == 0) {
      std::printf("Synchronization time: %.4f s\n", outputs.GetSyncTime());
    }
    outputs.CopyNewStateToCurrentState(inputs);
    ++cycles;
  }

  std::printf("Reached target after %d cycles (%.3f s)\n", cycles,
              cycles * kCycleTime);
  for (const auto& dof : inputs.GetDOFs()) {
    std::printf("  dof %d: pos=%.6f vel=%.6f acc=%.6f\n", dof.index,
                dof.position, dof.velocity, dof.acceleration);
  }
  return 0;
}
