#pragma once
#include <cg/types/amp.h>
#include <cg/types/box.h>

namespace cg::qa {
class prepare_nh {
public:
  nitro::vector<vec3r> const *r;
  nitro::vector<vec3r> *n, *h;
  box<real> const *box;
  nitro::vector<int> const *prev, *next;
  int num_particles;

public:
  void iter(int idx) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::qa