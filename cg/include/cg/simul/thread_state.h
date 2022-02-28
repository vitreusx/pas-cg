#pragma once
#include "kernels.h"
#include <cg/simul/dynamics.h>
#include <cg/utils/random.h>

namespace cg::simul {
class thread_state {
public:
  dynamics dyn;
  rand_gen gen;
  kernels ker;

public:
  thread_state() = default;
  explicit thread_state(int num_residues, uint64_t seed, kernels ker);

private:
  void fill_kernels();
};
} // namespace cg::simul