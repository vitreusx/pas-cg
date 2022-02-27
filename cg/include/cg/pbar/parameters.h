#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::pbar {
struct parameters {
  bool enabled;
  size_t width;
  quantity update_period;
  void link(ioxx::xyaml::proxy &proxy);
};
} // namespace cg::pbar