#pragma once
#include "data.h"

namespace cg::nl {
class verify {
public:
  nitro::vector<vec3r> const *r;
  box<real> const *box;
  bool *invalid, *first_time;
  int num_particles;
  nl_data const *data;

public:
  void operator()() const;
};
} // namespace cg::nl