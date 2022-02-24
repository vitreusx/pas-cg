#pragma once
#include <cg/simul/dynamics.h>
#include <cg/utils/random.h>

namespace cg::simul {
class thread_state {
public:
  dynamics dyn;
  rand_gen gen;

public:
  thread_state() = default;
  explicit thread_state(int num_residues, uint64_t seed);
};
} // namespace cg::simul