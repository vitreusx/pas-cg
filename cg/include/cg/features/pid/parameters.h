#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/xyaml.h>

namespace cg::pid {
struct parameters {
  bool enabled;

  enum class lambda_version_t { UNKNOWN, COSINE, ALGEBRAIC };
  lambda_version_t lambda_version;

  struct bb_t {
    quantity alpha, psi_0, r_min, depth;
    void connect(ioxx::xyaml_proxy &proxy);
  };
  bb_t bb_plus, bb_minus;

  struct ss_t {
    quantity alpha, psi_0;
    void connect(ioxx::xyaml_proxy &proxy);
  };
  ss_t ss;

  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg::pid