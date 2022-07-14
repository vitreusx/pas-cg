#pragma once
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg::dh {
struct const_dh_parameters {
  std::optional<bool> enabled;
  quantity permittivity;
  void load(ioxx::xyaml::node const &node);
};

struct rel_dh_parameters {
  std::optional<bool> enabled;
  quantity perm_factor;
  void load(ioxx::xyaml::node const &node);
};

struct parameters {
  bool enabled;
  std::string variant;
  quantity screening_dist;
  const_dh_parameters const_dh;
  rel_dh_parameters rel_dh;
  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::dh