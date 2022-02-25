#pragma once
#include <cg/types/amp.h>

namespace cg::simul {
class dynamics {
public:
  nitro::vector<vec3r> F;
  real V = 0.0;

public:
  dynamics() = default;
  explicit dynamics(int num_residues);
  void reset();
  void omp_reduce(dynamics &target);
  void sanity_check();
};
} // namespace cg::simul