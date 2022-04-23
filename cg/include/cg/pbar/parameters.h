#pragma once
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg::pbar {
struct parameters {
  bool enabled;
  size_t width;
  quantity update_period;
  void link(ioxx::xyaml::proxy &proxy);
};
} // namespace cg::pbar