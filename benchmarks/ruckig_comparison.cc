// Google Benchmark comparing Reflexxes Type IV against Ruckig.
//
// Both libraries receive identical, randomly sampled position-control problems
// (initial state, target state, kinematic limits). Two scenarios are measured:
//
//  - Calculate: every call receives a new problem, forcing a full trajectory
//    computation plus sampling of the first cycle.
//  - Step: the trajectory has already been computed and every call only
//    advances it by one control cycle (the steady-state online cost).
//
// Reflexxes Type IV has no target acceleration, so Ruckig targets are sampled
// with zero acceleration. Both use phase synchronization when possible.

#include <benchmark/benchmark.h>

#include <array>
#include <cstddef>
#include <random>
#include <vector>

#include "intrinsic/icon/reflexxes/reflexxes_api.h"
#include "ruckig/ruckig.hpp"

namespace rx = intrinsic::reflexxes;

namespace {

constexpr double kCycleTime = 0.001;
constexpr std::size_t kNumProblems = 1024;

struct DofProblem {
  double position, velocity, acceleration;
  double target_position, target_velocity;
  double max_velocity, max_acceleration, max_jerk;
};

template <std::size_t N>
using Problem = std::array<DofProblem, N>;

template <std::size_t N>
std::vector<Problem<N>> MakeProblems() {
  std::mt19937_64 gen(42);
  auto uniform = [&](double lo, double hi) {
    return std::uniform_real_distribution<double>(lo, hi)(gen);
  };
  std::vector<Problem<N>> problems(kNumProblems);
  for (auto& problem : problems) {
    for (auto& dof : problem) {
      dof.max_velocity = uniform(1.0, 3.0);
      dof.max_acceleration = uniform(1.0, 5.0);
      dof.max_jerk = uniform(5.0, 20.0);
      dof.position = uniform(-4.0, 4.0);
      dof.velocity = uniform(-0.5, 0.5) * dof.max_velocity;
      dof.acceleration = uniform(-0.5, 0.5) * dof.max_acceleration;
      dof.target_position = uniform(-4.0, 4.0);
      dof.target_velocity = uniform(-0.5, 0.5) * dof.max_velocity;
    }
  }
  return problems;
}

// --- Reflexxes --------------------------------------------------------------

template <std::size_t N>
void SetReflexxesInputs(const Problem<N>& problem, rx::PositionInputs& inputs) {
  for (std::size_t i = 0; i < N; ++i) {
    const DofProblem& p = problem[i];
    rx::Inputs::DOF& dof = inputs.GetDOFs()[i];
    dof.index = static_cast<int>(i);
    dof.position = p.position;
    dof.velocity = p.velocity;
    dof.acceleration = p.acceleration;
    dof.target_position = p.target_position;
    dof.target_velocity = p.target_velocity;
    dof.max_position = 100.0;
    dof.min_position = -100.0;
    dof.max_velocity = p.max_velocity;
    dof.min_velocity = -p.max_velocity;
    dof.max_acceleration = p.max_acceleration;
    dof.min_acceleration = -p.max_acceleration;
    dof.max_jerk = p.max_jerk;
    dof.min_jerk = -p.max_jerk;
  }
}

bool IsReflexxesOk(rx::Status status) {
  return status == rx::Status::kWorking ||
         status == rx::Status::kFinalStateReached;
}

template <std::size_t N>
void BM_Reflexxes_Calculate(benchmark::State& bm) {
  const auto problems = MakeProblems<N>();
  rx::State state(N, kCycleTime);
  rx::PositionInputs inputs(N, kCycleTime);
  rx::PositionOutputs outputs(N, kCycleTime);
  const rx::PositionFlags flags;

  std::size_t idx = 0;
  double total_duration = 0.0;
  int64_t failures = 0;
  for (auto _ : bm) {
    SetReflexxesInputs(problems[idx], inputs);
    const rx::Status status = rx::ComputePosition(inputs, flags, outputs, state);
    benchmark::DoNotOptimize(outputs);
    failures += !IsReflexxesOk(status);
    total_duration += outputs.GetSyncTime();
    idx = (idx + 1) % kNumProblems;
  }
  bm.counters["avg_duration_s"] =
      benchmark::Counter(total_duration, benchmark::Counter::kAvgIterations);
  bm.counters["failures"] = static_cast<double>(failures);
}

template <std::size_t N>
void BM_Reflexxes_Step(benchmark::State& bm) {
  const auto problems = MakeProblems<N>();
  rx::State state(N, kCycleTime);
  rx::PositionInputs inputs(N, kCycleTime);
  rx::PositionOutputs outputs(N, kCycleTime);
  const rx::PositionFlags flags;

  std::size_t idx = 0;
  auto reset = [&] {
    SetReflexxesInputs(problems[idx], inputs);
    idx = (idx + 1) % kNumProblems;
    rx::ComputePosition(inputs, flags, outputs, state);
    outputs.CopyNewStateToCurrentState(inputs);
  };
  reset();

  int64_t failures = 0;
  for (auto _ : bm) {
    const rx::Status status = rx::ComputePosition(inputs, flags, outputs, state);
    benchmark::DoNotOptimize(outputs);
    failures += !IsReflexxesOk(status);
    outputs.CopyNewStateToCurrentState(inputs);
    if (status != rx::Status::kWorking) {
      bm.PauseTiming();
      reset();
      bm.ResumeTiming();
    }
  }
  bm.counters["failures"] = static_cast<double>(failures);
}

// --- Ruckig -----------------------------------------------------------------

template <std::size_t N>
void SetRuckigInput(const Problem<N>& problem,
                    ruckig::InputParameter<N>& input) {
  input.synchronization = ruckig::Synchronization::Phase;
  for (std::size_t i = 0; i < N; ++i) {
    const DofProblem& p = problem[i];
    input.current_position[i] = p.position;
    input.current_velocity[i] = p.velocity;
    input.current_acceleration[i] = p.acceleration;
    input.target_position[i] = p.target_position;
    input.target_velocity[i] = p.target_velocity;
    input.target_acceleration[i] = 0.0;
    input.max_velocity[i] = p.max_velocity;
    input.max_acceleration[i] = p.max_acceleration;
    input.max_jerk[i] = p.max_jerk;
  }
}

bool IsRuckigOk(ruckig::Result result) {
  return result == ruckig::Result::Working ||
         result == ruckig::Result::Finished;
}

template <std::size_t N>
void BM_Ruckig_Calculate(benchmark::State& bm) {
  const auto problems = MakeProblems<N>();
  ruckig::Ruckig<N> otg(kCycleTime);
  ruckig::InputParameter<N> input;
  ruckig::OutputParameter<N> output;

  std::size_t idx = 0;
  double total_duration = 0.0;
  int64_t failures = 0;
  for (auto _ : bm) {
    SetRuckigInput(problems[idx], input);
    const ruckig::Result result = otg.update(input, output);
    benchmark::DoNotOptimize(output);
    failures += !IsRuckigOk(result);
    total_duration += output.trajectory.get_duration();
    idx = (idx + 1) % kNumProblems;
  }
  bm.counters["avg_duration_s"] =
      benchmark::Counter(total_duration, benchmark::Counter::kAvgIterations);
  bm.counters["failures"] = static_cast<double>(failures);
}

template <std::size_t N>
void BM_Ruckig_Step(benchmark::State& bm) {
  const auto problems = MakeProblems<N>();
  ruckig::Ruckig<N> otg(kCycleTime);
  ruckig::InputParameter<N> input;
  ruckig::OutputParameter<N> output;

  std::size_t idx = 0;
  auto reset = [&] {
    SetRuckigInput(problems[idx], input);
    idx = (idx + 1) % kNumProblems;
    otg.update(input, output);
    output.pass_to_input(input);
  };
  reset();

  int64_t failures = 0;
  for (auto _ : bm) {
    const ruckig::Result result = otg.update(input, output);
    benchmark::DoNotOptimize(output);
    failures += !IsRuckigOk(result);
    output.pass_to_input(input);
    if (result != ruckig::Result::Working) {
      bm.PauseTiming();
      reset();
      bm.ResumeTiming();
    }
  }
  bm.counters["failures"] = static_cast<double>(failures);
}

}  // namespace

#define REGISTER_DOFS(N)                              \
  BENCHMARK(BM_Reflexxes_Calculate<N>)->Name("Calculate/Reflexxes/" #N "dof"); \
  BENCHMARK(BM_Ruckig_Calculate<N>)->Name("Calculate/Ruckig/" #N "dof");       \
  BENCHMARK(BM_Reflexxes_Step<N>)->Name("Step/Reflexxes/" #N "dof");           \
  BENCHMARK(BM_Ruckig_Step<N>)->Name("Step/Ruckig/" #N "dof");

REGISTER_DOFS(1)
REGISTER_DOFS(3)
REGISTER_DOFS(6)
REGISTER_DOFS(7)

BENCHMARK_MAIN();
