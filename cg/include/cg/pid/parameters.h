#pragma once
#include "lambda.h"
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::pid {
struct parameters {
  bool enabled;

  pid::lambda_version lambda_version;

  struct bb_t {
    quantity alpha, psi_0, r_min, depth;
    void load(ioxx::xyaml::node const &node);
  };
  bb_t bb_plus, bb_minus;

  struct ss_t {
    quantity alpha, psi_0;
    void load(ioxx::xyaml::node const &node);
  };
  ss_t ss;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::pid