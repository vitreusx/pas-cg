#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/xyaml.h>

namespace cg::pbar {
struct parameters {
  bool enabled;
  size_t width;
  quantity update_period;
  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg::pbar