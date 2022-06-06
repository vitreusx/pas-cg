#pragma once
#include "pair.h"
#include <cg/types/box.h>
#include <optional>

namespace cg::nl {
class data {
public:
  vect::vector<pair> native, non_native;
  real orig_pad;
  std::optional<real> fixed_cutoff;
  box<real> orig_box;
  vect::vector<vec3r> orig_r;

  inline bool in_range(real dist, real cutoff) const {
    if (fixed_cutoff.has_value())
      cutoff = fixed_cutoff.value();
    return dist < cutoff + orig_pad;
  }
};
} // namespace cg::nl