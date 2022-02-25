#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/xyaml.h>

namespace cg::gen {
struct parameters {
  quantity total_time;
  uint64_t seed, num_of_threads;
  bool debug_mode;

  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg::gen