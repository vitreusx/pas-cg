#pragma once
#include "contact_type.h"
#include <cg/base_forces/sink_lj.h>
#include <nitro/nitro.h>

namespace cg::qa {
class lj_variants {
public:
  inline decltype(auto) operator[](contact_type const &type) const {
    return variants[(int16_t)type];
  }

private:
  nitro::vector<sink_lj> variants;
};
} // namespace cg::qa