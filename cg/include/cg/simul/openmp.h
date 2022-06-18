#pragma once
#include "parameters.h"
#include "state.h"

namespace cg::simul {
class openmp {
public:
  explicit openmp(state const &st);
  void main();

private:
  state st;
};
} // namespace cg::simul