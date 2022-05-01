#pragma once
#include "../at.h"
#include "../const_at.h"
#include "../const_view.h"
#include "../expr.h"
#include "../lane_at.h"
#include "../lane_const_at.h"
#include "../view.h"
#include "decl.h"
#include <cstddef>

namespace nitro {
template <typename Types, typename Allocs, typename Idx, size_t... ISeq>
class par_vector {
public:
  template <size_t I>
  using slice = vector<typename Types::template ith_type<I>,
                       typename Allocs::template ith_type<I>, Idx>;

  par_vector() : par_vector(Allocs()){};

  explicit par_vector(Allocs const &allocs)
      : slices{slice<ISeq>(allocs.template get<ISeq>())...} {};

  explicit par_vector(Idx n, Types const &init = Types(),
                      Allocs allocs = Allocs())
      : slices{slice<ISeq>(n, init.template get<ISeq>(),
                           allocs.template get<ISeq>())...} {}

  template <typename E>
  explicit par_vector(Idx n, ind_expr<E> const &e, Allocs allocs = Allocs())
      : slices{slice<ISeq>(n, e.template get<ISeq>(),
                           allocs.template get<ISeq>())...} {}

  Idx size() const {
    return slices.template get<0>().size();
  }
  Idx capacity() const {
    return slices.template get<0>().capacity();
  }

  at_expr<Types> operator[](Idx idx) {
    return get_view()[idx];
  }

  const_at_expr<Types> operator[](Idx idx) const {
    return get_view()[idx];
  }

  at_expr<Types> at(Idx idx) {
    return get_view().at(idx);
  }

  const_at_expr<Types> at(Idx idx) const {
    return get_view().at(idx);
  }

  template <size_t N> def_lane_at<Types, N> lane_at(Idx idx) {
    return get_view().lane_at(idx);
  }

  template <size_t N> lane<Types, N> lane_at(Idx idx) const {
    return get_view().lane_at(idx);
  }

  view<Types, Idx> get_view() {
    return par_view<Types, Idx, ISeq...>(
        slices.template get<ISeq>()...);
  }

  operator view<Types, Idx>() {
    return get_view();
  }

  const_view<Types, Idx> get_view() const {
    return par_const_view<Types, Idx, ISeq...>(
        slices.template get<ISeq>()...);
  }

  operator const_view<Types, Idx>() const {
    return get_view();
  }

  void clear() {
    (..., slices.template get<ISeq>().clear());
  }

  void reserve(Idx new_capacity) {
    (..., slices.template get<ISeq>().reserve(new_capacity));
  }

  void resize(Idx new_size) {
    resize(new_size, Types());
  }

  template <typename E> void resize(Idx new_size, ind_expr<E> const &e) {
    (..., slices.template get<ISeq>().resize(new_size, e.template get<ISeq>()));
  }

  void shrink(Idx new_size) {
    (..., slices.template get<ISeq>().shrink(new_size));
  }

  template <typename E> void push_back(ind_expr<E> const &e) {
    (..., slices.template get<ISeq>().push_back(e.template get<ISeq>()));
  }

  template <typename... Args> at_expr<Types> emplace_back(Args &&...args) {
    push_back(Types(args...));
    return at(size() - 1);
  }

  iterator<Types> begin() {
    return get_view().begin();
  }

  const_iterator<Types> begin() const {
    return get_view().begin();
  }

  iterator<Types> end() {
    return get_view().end();
  }

  const_iterator<Types> end() const {
    return get_view().end();
  }

private:
  tuple<slice<ISeq>...> slices;
};
} // namespace nitro
