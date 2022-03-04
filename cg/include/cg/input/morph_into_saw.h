#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::input {
struct morph_into_saw_t {
  bool perform;
  std::optional<quantity> bond_distance;
  quantity residue_density, intersection_at;
  bool infer_box;
  int num_of_retries;
  void link(ioxx::xyaml::proxy &proxy);
};
} // namespace cg::input