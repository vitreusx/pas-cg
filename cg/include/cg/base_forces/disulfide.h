#pragma once
#include "harmonic.h"
#include "lj.h"
#include <cg/files/files.h>
#include <variant>

namespace cg {
class disulfide_force {
public:
  std::variant<std::monostate, harmonic, lj> force;

  __host__ __device__ auto operator()(real norm_inv) const {
    if (std::holds_alternative<harmonic>(force)) {
      return std::get<harmonic>(force)(1.0 / norm_inv);
    } else if (std::holds_alternative<lj>(force)) {
      return std::get<lj>(force)(norm_inv);
    } else {
      return vect::tuple<real, real>((real)0, (real)0);
    }
  }
};
} // namespace cg