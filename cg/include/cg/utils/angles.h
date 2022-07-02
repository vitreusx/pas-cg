#pragma once
#include "math.h"

namespace cg {
template <typename V1, typename V2, typename V3>
auto bond_angle(V1 const &r1, V2 const &r2, V3 const &r3) {
  auto u1 = r2 - r1, u2 = r2 - r3;
  return acos(dot(u1, u2) / (norm(u1) * norm(u2)));
}

template <typename V1, typename V2, typename V3, typename V4>
auto dihedral_angle(V1 const &r1, V2 const &r2, V3 const &r3, V4 const &r4) {
  auto u1 = r2 - r1, u2 = r3 - r2, u3 = r4 - r3;
  auto v1 = cross(u1, u2), v2 = cross(u2, u3);
  auto phi = acos(dot(v1, v2) / (norm(v1) * norm(v2)));
  if (dot(v1, u3) < 0)
    phi = -phi;
  return phi;
}
} // namespace cg