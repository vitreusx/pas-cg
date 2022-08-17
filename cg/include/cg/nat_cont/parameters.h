#pragma once
#include <cg/base_forces/spec.h>
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg::nat_cont {
struct parameters {
  bool enabled;
  quantity lj_depth;
  std::optional<force_spec> ss_force;

  struct unfolding_study_t {
    bool early_stopping, measure_times;
    void load(ioxx::xyaml::node const &n);
  };
  unfolding_study_t unfolding_study;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::nat_cont