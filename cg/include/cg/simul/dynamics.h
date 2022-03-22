#pragma once
#include <cg/types/amp.h>

namespace cg::simul {
class thread;

class dynamics {
public:
  nitro::vector<vec3r> F;
  real V = 0.0;

public:
  dynamics() = default;
  explicit dynamics(int num_residues);

  void reset();
  void omp_reset();

  void reduce(dynamics &target);
  void omp_reduce(dynamics &target);
  void omp_reduce_v2(dynamics &target, thread const &thr);
};
} // namespace cg::simul