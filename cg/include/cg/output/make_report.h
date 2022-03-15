#pragma once
#include "hook.h"
#include <cg/types/amp.h>
#include <vector>

namespace cg::out {

class make_report {
public:
  real stats_period, file_period;
  real *t;
  report_data *state;
  std::vector<hook const *> const *hooks;

public:
  void operator()() const;
};
} // namespace cg::out