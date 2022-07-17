#pragma once
#include <cg/amino/amino_acid.h>
#include <cg/files/files.h>
#include <cg/types/amp.h>
#include <cg/utils/hash.h>
#include <cg/utils/quantity.h>

namespace cg {
class shifted_lj {
public:
  real depth, r_min;

  shifted_lj() = default;

  inline shifted_lj(real depth, real r_min) : depth{depth}, r_min{r_min} {};

  inline auto operator()(real r_inv) const {
    auto s = r_inv * r_min, s6 = ipow<6>(s), s12 = s6 * s6;
    auto V = depth * (s12 - (real)2.0 * s6 + (real)1.0);
    auto dV_dr = (real)12.0 * depth * r_inv * (s6 - s12);
    return std::make_tuple(V, dV_dr);
  }

private:
};

struct shifted_lj_specs {
  std::optional<quantity> depth, r_min;
  operator shifted_lj() const;
  void load(ioxx::xyaml::node const &node);
};

struct ss_shifted_lj_specs {
  std::unordered_map<std::pair<amino_acid, amino_acid>, shifted_lj_specs>
      ss_specs;
  void load(ioxx::xyaml::node const &node);
};
} // namespace cg