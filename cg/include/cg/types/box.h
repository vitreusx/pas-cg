#pragma once
#include <cg/types/vec3.h>
#include <cg/utils/math.h>

namespace cg {
template <typename U> class box {
public:
  vec3<U> cell, cell_inv;

  inline box() : cell{vec3<U>::Zero()}, cell_inv{vec3<U>::Zero()} {};

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

  template <typename E> inline auto wrap(E const &e) const {
    vec3<U> in_box_ = e;
    in_box_.x() -= cell.x() * nearbyint(e.x() * cell_inv.x());
    in_box_.y() -= cell.y() * nearbyint(e.y() * cell_inv.y());
    in_box_.z() -= cell.z() * nearbyint(e.z() * cell_inv.z());
    return in_box_;
  }

  template <typename E, typename F>
  inline auto wrap(E const &e, F const &f) const {
    vec3<U> diff = f - e;
    return wrap(diff);
  }
};
} // namespace cg
