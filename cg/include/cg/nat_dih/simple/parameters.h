#pragma once
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg::snd {
struct parameters {
  bool enabled;
  quantity CDH;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::snd