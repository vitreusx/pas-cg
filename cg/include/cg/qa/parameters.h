#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::qa {
struct parameters {
  bool enabled;
  quantity phase_dur;
  double breaking_factor, min_cos_hr, min_cos_hh, max_cos_nr;
  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::qa