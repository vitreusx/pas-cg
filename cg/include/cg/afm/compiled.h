#pragma once
#include "force/tip.h"
#include "parameters.h"
#include "vel/tip.h"

namespace cg::afm {
class compiled_tips {
public:
  nitro::vector<force::tip> force;
  nitro::vector<vel::tip> vel;
  nitro::vector<int> pulled_chains;
};
} // namespace cg::afm