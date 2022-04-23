#pragma once
#include "../const_at.h"
#include "../lane_const_at.h"
#include "../const_iterator.h"
#include "decl.h"

namespace nitro {
template <typename T, typename Idx> class def_const_view {
public:
  def_const_view() = default;
  explicit def_const_view(T const *data, Idx n) : data{data}, n{n} {};

  const_at_expr<T> operator[](Idx idx) const { return data[idx]; }

  const_at_expr<T> at(Idx idx) const { return (*this)[idx]; }

  template <size_t N> lane_const_at_expr<T, N> lane_at(Idx idx) const {
    return lane_const_at_expr<T, N>(data[N * idx]);
  }

  template <size_t N> Idx num_lanes() const { return size() / N; }

  Idx size() const { return n; }

  const_iterator<T> begin() const {
    return def_const_iterator<T>(data);
  }

  const_iterator<T> end() const {
    return def_const_iterator<T>(data + n);
  }

protected:
  T const *data;
  Idx n;
};
} // namespace nitro