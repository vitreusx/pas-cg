#pragma once
#include <cg/files/files.h>
#include <cg/types/vec3.h>
#include <cg/utils/quantity.h>
#include <variant>
#include <vector>

namespace cg::afm {
struct parameters {
  bool enabled;
  quantity H1, H2;

  enum class tip_type { CONST_VEL, CONST_FORCE };

  struct single_res_t {
    int res_idx;
    tip_type type;
    vec3<quantity> dir;
  };

  struct pulled_apart_t {
    int chain_idx;
    tip_type type;
    quantity mag;
  };

  using tip_t = std::variant<single_res_t, pulled_apart_t>;
  std::vector<tip_t> tips;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::afm