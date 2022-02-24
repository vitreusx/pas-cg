#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/xyaml.h>

namespace cg::qa {
struct parameters {
  bool enabled;
  quantity phase_dur;
  double breaking_factor, min_cos_hr, min_cos_hh, max_cos_nr;
  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg::qa