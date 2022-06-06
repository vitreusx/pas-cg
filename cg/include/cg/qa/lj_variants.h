#pragma once
#include "contact_type.h"
#include <cg/base_forces/lj_variants.h>
#include <cg/base_forces/sink_lj.h>
#include <cg/vect/vect.h>

namespace cg::qa {
class lj_variants {
public:
  lj_variants() = default;
  explicit lj_variants(cg::lj_variants const &variants);

  inline decltype(auto) operator[](contact_type const &type) const {
    return variants[(int16_t)type];
  }

private:
  vect::vector<sink_lj> variants;
};
} // namespace cg::qa