#pragma once
#include <cg/base_forces/disulfide.h>
#include <cg/base_forces/lj.h>
#include <cg/base_forces/sink_lj.h>
#include <cg/base_forces/spec.h>
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg::qa {

struct parameters {
  bool enabled, include4;
  quantity phase_dur;
  double breaking_factor, min_cos_hr, min_cos_hh, max_cos_nr,
      formation_tolerance;

  std::string def_force_variant;
  force_spec bb, bs;
  ss_force_spec ss;

  struct ss_spec_crit_t {
    bool enabled;
    int max_neigh_count;
    quantity def_dist, max_dist_dev, neigh_radius;
    void load(ioxx::xyaml::node const &node);
  };

  struct disulfide_t {
    force_spec force;
    ss_spec_crit_t spec_crit;
    void load(ioxx::xyaml::node const &node);
  };
  std::optional<disulfide_t> disulfide;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::qa