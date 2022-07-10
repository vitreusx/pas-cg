#pragma once
#include "lambda.h"
#include <cg/amino/amino_acid.h>
#include <cg/base_forces/lj.h>
#include <cg/base_forces/sink_lj.h>
#include <cg/base_forces/spec.h>
#include <cg/files/files.h>
#include <cg/utils/quantity.h>
#include <optional>
#include <unordered_map>
#include <variant>

namespace cg::pid {
struct parameters {
  bool enabled, include4;

  struct lambda_params {
    quantity alpha, psi_0;
    void load(ioxx::xyaml::node const &node);
  };

  pid::lambda_version lambda_variant;
  lambda_params bb_plus_lambda, bb_minus_lambda, ss_lambda;

  std::string def_force_variant;
  force_spec bb_plus_force, bb_minus_force;
  ss_force_spec ss_force;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::pid