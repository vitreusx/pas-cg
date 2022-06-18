#pragma once
#include "pair.h"
#include <cg/types/box.h>
#include <optional>

namespace cg::nl {
class data {
public:
  vect::vector<pair> native, non_native;
  real orig_pad;
  box<real> orig_box;
  vect::vector<vec3r> orig_r;
};
} // namespace cg::nl