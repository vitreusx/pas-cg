#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/xyaml.h>

namespace cg::snd {
struct parameters {
  bool enabled;
  quantity CDH;
  
  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg::snd