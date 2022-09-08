#pragma once
#include "pair.h"
#include <cg/sbox/pbc.h>
#include <optional>

namespace cg::nl {
class data {
public:
  vect::vector<pair> pairs;
  real orig_pad;
  sbox::pbc<real> orig_pbc;
  vect::vector<vec3r> orig_r;
};
} // namespace cg::nl