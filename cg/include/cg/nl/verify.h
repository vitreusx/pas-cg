#pragma once
#include "data.h"

namespace cg::nl {
class verify {
public:
  nitro::const_view<vec3r> r;
  box<real> const *simul_box;
  bool *invalid, *first_time;
  int num_particles;
  data const *nl_data;

public:
  void operator()() const;
};
} // namespace cg::nl