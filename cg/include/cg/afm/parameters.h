#pragma once
#include <cg/files/files.h>
#include <cg/types/vec3.h>
#include <cg/utils/quantity.h>
#include <variant>
#include <vector>

namespace cg::afm {
struct parameters {
  bool perform;
  std::string type;

  struct pull_rel_t {
    quantity time;
    void load(ioxx::xyaml::node const& n);
  };
  pull_rel_t pull_rel;

  struct fafm_t {
    quantity force;
    void load(ioxx::xyaml::node const &n);
  };

  struct vafm_t {
    quantity vel, H1, H2;
    void load(ioxx::xyaml::node const &n);
  };

  struct tip_params_t {
    std::string type;
    fafm_t force_afm;
    vafm_t vel_afm;
    void load(ioxx::xyaml::node const &n);
  };
  tip_params_t tip_params;

  void load(ioxx::xyaml::node const &n);
};
} // namespace cg::afm