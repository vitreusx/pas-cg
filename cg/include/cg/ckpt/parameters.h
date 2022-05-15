#pragma once
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg::ckpt {
struct parameters {
  bool enabled;
  std::string path_fmt;
  quantity every;

  void load(ioxx::xyaml::node const& node);
};
}