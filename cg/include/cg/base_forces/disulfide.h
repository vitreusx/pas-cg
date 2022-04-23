#pragma once
#include "harmonic.h"
#include "lj.h"
#include <cg/files/files.h>
#include <variant>

namespace cg {
class disulfide_force {
public:
  std::variant<std::monostate, harmonic, lj> force;

  std::tuple<real, real> operator()(real norm_inv) const {
    if (std::holds_alternative<harmonic>(force)) {
      return std::get<harmonic>(force)(1.0 / norm_inv);
    } else if (std::holds_alternative<lj>(force)) {
      return std::get<lj>(force)(norm_inv);
    } else {
      return std::make_tuple((real)0, (real)0);
    }
  }

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg