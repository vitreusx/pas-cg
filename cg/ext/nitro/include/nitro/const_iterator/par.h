#pragma once
#include "../const_at.h"
#include "../tuple.h"
#include "def.h"
#include <cstddef>

namespace nitro {
template <typename Types, size_t... ISeq> class par_const_iterator {
public:
  template <size_t I>
  using slice = const_iterator<typename Types::template ith_type<I>>;

  explicit par_const_iterator(slice<ISeq> const &...slices)
      : slices{slices...} {};

  const_at_expr<Types> operator*() const {
    return {*slices.template get<ISeq>()...};
  }

  par_const_iterator &operator++() {
    (..., ++slices.template get<ISeq>());
    return *this;
  }

  par_const_iterator operator++(int) {
    return par_const_iterator(slices.template get<ISeq>()++...);
  }

  template <typename Idx>
  par_const_iterator operator+(Idx const &offset) const {
    return par_const_iterator(slices.template get<ISeq>() + offset...);
  }

  template <typename Idx> par_const_iterator &operator+=(Idx const &offset) {
    (..., advance<ISeq>(offset));
    return *this;
  }

  par_const_iterator &operator--() {
    (..., --slices.template get<ISeq>());
    return *this;
  }

  par_const_iterator operator--(int) {
    return par_const_iterator(slices.template get<ISeq>()--...);
  }

  template <typename Idx>
  par_const_iterator operator-(Idx const &offset) const {
    return par_const_iterator(slices.template get<ISeq>() - offset...);
  }

  template <typename Idx> par_const_iterator &operator-=(Idx const &offset) {
    (..., backtrack<ISeq>(offset));
    return *this;
  }

  bool operator<(par_const_iterator const &other) const {
    return slices.template get<0>() < other.slices.template get<0>();
  }

  bool operator<=(par_const_iterator const &other) const {
    return slices.template get<0>() <= other.slices.template get<0>();
  }

  bool operator>(par_const_iterator const &other) const {
    return slices.template get<0>() > other.slices.template get<0>();
  }

  bool operator>=(par_const_iterator const &other) const {
    return slices.template get<0>() >= other.slices.template get<0>();
  }

  bool operator==(par_const_iterator const &other) const {
    return slices.template get<0>() == other.slices.template get<0>();
  }

  bool operator!=(par_const_iterator const &other) const {
    return slices.template get<0>() != other.slices.template get<0>();
  }

private:
  template <size_t I, typename Idx> void advance(Idx const &offset) {
    slices.template get<I>() += offset;
  }

  template <size_t I, typename Idx> void backtrack(Idx const &offset) {
    slices.template get<I>() -= offset;
  }

  tuple<slice<ISeq>...> slices;
};
} // namespace nitro