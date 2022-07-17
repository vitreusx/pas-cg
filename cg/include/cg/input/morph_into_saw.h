#pragma once
#include <cg/files/files.h>
#include <cg/utils/quantity.h>
#include <variant>

namespace cg::input {
struct morph_into_saw_t {
  struct start_box_params {
    std::optional<quantity> density;
    std::optional<quantity> size;
  };
  struct start_box_origin {};
  std::variant<std::monostate, start_box_params, start_box_origin> start_box;

  bool perform;
  std::optional<quantity> bond_distance;
  quantity intersection_at;
  bool with_pbc;
  int num_of_retries;
  void link(ioxx::xyaml::proxy &proxy);
};
} // namespace cg::input