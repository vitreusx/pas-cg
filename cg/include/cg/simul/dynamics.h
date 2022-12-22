#pragma once
#include <cg/types/amp.h>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <stack>
#include <utility>

namespace cg::simul {
class thread;
class thread_team;

class dynamics {
public:
  vect::vector<vec3r> F, solid_wall_F, harmonic_wall_F, lj_wall_F;
  real V = 0.0;

public:
  dynamics() = default;
  explicit dynamics(int num_residues);

  void reset();
  void omp_reset();

  void reduce(dynamics &target);
  void omp_reduce(dynamics &target);
  void omp_reduce_v2(dynamics &target, thread const &thr);
  void omp_reduce_v3(dynamics &target, thread_team const &team);

  struct v4_shared {
    int num_threads;

    struct section {
      section() = default;

      int beg = 0, end = 0;
      std::shared_ptr<std::mutex> mtx;
    };
    std::vector<section> sections;
  };
  struct v4_priv {
    std::vector<int> sec_order;
  };
  void omp_reduce_v4(dynamics &target, v4_shared &shared, v4_priv &priv);
};
} // namespace cg::simul