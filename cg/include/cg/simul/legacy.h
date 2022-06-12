#pragma once
#include "kernels.h"
#include "state.h"

namespace cg::simul {
class legacy {
public:
  explicit legacy(state &st);

  void main();

public:
  state &st;
};
} // namespace cg::simul