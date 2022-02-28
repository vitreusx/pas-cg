#pragma once
#include "hook.h"
#include <cg/types/amp.h>
#include <vector>

namespace cg::out {

class make_report {
public:
  real period;
  real *t, *last_t;
  report_state *state;
  std::vector<hook const *> const *hooks;

public:
  void operator()() const;
};
} // namespace cg::out