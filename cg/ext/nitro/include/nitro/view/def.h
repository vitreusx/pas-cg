#pragma once
#include "../at.h"
#include "../const_iterator.h"
#include "../const_view.h"
#include "../iterator.h"
#include "../lane_at.h"
#include "nitro/at/decl.h"

namespace nitro {
template <typename T, typename Idx> class def_view {
public:
  def_view() = default;
  explicit def_view(T *data, Idx n) : data{data}, n{n} {};

  at_expr<T> operator[](Idx idx) const { return data[idx]; }

  at_expr<T> at(Idx idx) const { return (*this)[idx]; }

  template <size_t N> lane_at_expr<T, N> lane_at(Idx idx) const {
    return lane_at<T, N>(data[N * idx]);
  }

  template <size_t N> Idx num_lanes() const { return size() / N; }

  Idx size() const { return n; }

  operator const_view<T, Idx>() const { return const_view<T, Idx>(data, n); }

  iterator<T> begin() { return def_iterator<T>(data); }

  const_iterator<T> begin() const { return def_const_iterator<T>(data); }

  iterator<T> end() { return def_iterator<T>(data + n); }

  const_iterator<T> end() const { return def_const_iterator<T>(data + n); }

private:
  T *data;
  Idx n;
};
} // namespace nitro