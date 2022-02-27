#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::gen {
struct parameters {
  quantity total_time;
  uint64_t seed, num_of_threads;
  bool debug_mode;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::gen