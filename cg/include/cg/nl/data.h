#pragma once
#include "pair.h"
#include <cg/types/box.h>

namespace cg::nl {
class data {
public:
  nitro::vector<pair> native, non_native;
  real orig_pad;
  box<real> orig_box;
  nitro::vector<vec3r> orig_r;
};
} // namespace cg::nl