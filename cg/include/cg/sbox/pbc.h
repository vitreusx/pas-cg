#pragma once
#include <cg/types/vec3.h>
#include <cg/utils/math.h>

namespace cg::sbox {
template <typename U>
class pbc {
public:
  vec3<U> cell, cell_inv;

  inline pbc() : cell{vec3<U>::Zero()}, cell_inv{vec3<U>::Zero()} {};

  inline void set_cell(vec3<U> const &new_cell) {
    this->cell = new_cell;

    cell_inv = vec3<U>::Zero();
    if (cell.x() != 0.0)
      cell_inv.x() = (U)1.0 / cell.x();
    if (cell.y() != 0.0)
      cell_inv.y() = (U)1.0 / cell.y();
    if (cell.z() != 0.0)
      cell_inv.z() = (U)1.0 / cell.z();
  }

  inline void set_cell_inv(vec3<U> const &new_cell_inv) {
    this->cell_inv = new_cell_inv;

    cell = vec3<U>::Zero();
    if (cell_inv.x() != 0.0)
      cell.x() = (U)1.0 / cell_inv.x();
    if (cell_inv.y() != 0.0)
      cell.y() = (U)1.0 / cell_inv.y();
    if (cell_inv.z() != 0.0)
      cell.z() = (U)1.0 / cell_inv.z();
  }

  template <typename V>
  __host__ __device__ inline auto wrap(V v) const {
    v.x() -= cell.x() * round(v.x() * cell_inv.x());
    v.y() -= cell.y() * round(v.y() * cell_inv.y());
    v.z() -= cell.z() * round(v.z() * cell_inv.z());
    return v;
  }

  template <typename V, typename E, typename F>
  __host__ __device__ inline auto wrap(E const &e, F const &f) const {
    auto diff = f - e;
    return wrap<V>(diff);
  }
};
} // namespace cg::sbox