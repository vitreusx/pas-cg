#pragma once
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg::input {
struct morph_into_saw_t {
  bool perform;
  std::optional<quantity> bond_distance;
  quantity init_box_density, intersection_at;
  bool with_pbc;
  int num_of_retries;
  void link(ioxx::xyaml::proxy &proxy);
};
} // namespace cg::input