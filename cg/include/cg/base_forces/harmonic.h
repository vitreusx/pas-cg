#pragma once
#include <cg/files/files.h>
#include <cg/types/amp.h>
#include <cg/utils/quantity.h>

namespace cg {
class harmonic {
public:
  real H1, H2;
  real nat_r;

  harmonic() = default;

  __host__ __device__ inline harmonic(real H1, real H2, real nat_r)
      : H1{H1}, H2{H2}, nat_r{nat_r} {};

public:
  __host__ __device__ inline auto operator()(real r) const {
    auto dr = r - nat_r, dr2 = dr * dr;
    auto V = (real)0.5 * dr2 * (H1 + H2 * dr2);
    auto dV_dr = dr * (H1 + (real)2.0 * H2 * dr2);
    dV_dr = clamp<real>(dV_dr, -1.0e3, 1.0e3);
    return vect::tuple<real, real>(V, dV_dr);
  }

  static inline auto cutoff(real nat_r) {
    return (real)2.0 * nat_r;
  }
};

struct harmonic_specs {
  std::optional<quantity> H1, H2, nat_r;
  void load(ioxx::xyaml::node const &node);
  operator harmonic() const;
};
} // namespace cg