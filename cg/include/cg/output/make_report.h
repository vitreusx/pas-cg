#pragma once
#include "hook.h"
#include <cg/types/amp.h>
#include <vector>

namespace cg::output {
class make_report : public hook {
public:
  real period, last_t;
  real *t;
  std::filesystem::path output_dir;
  int *ord;
  std::vector<hook const *> hooks;

public:
  void operator()() const;
};
} // namespace cg::output