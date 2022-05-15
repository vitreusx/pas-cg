#pragma once
#include <cg/simul/state.h>

namespace cg::ckpt {
class make_checkpoint {
public:
  std::string path_fmt;
  real *last_t, every;
  simul::state const *st;

public:
  void operator()() const;
};
} // namespace cg::ckpt