#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::snd {
struct parameters {
  bool enabled;
  quantity CDH;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::snd