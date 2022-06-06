#pragma once
#include "force/tip.h"
#include "parameters.h"
#include "vel/tip.h"

namespace cg::afm {
class compiled_tips {
public:
  vect::vector<force::tip> force;
  vect::vector<vel::tip> vel;
  vect::vector<int> pulled_chains;
};
} // namespace cg::afm