#pragma once
#include "../at/decl.h"
#include "../tuple.h"
#include "def.h"
#include <cstddef>

namespace nitro {
template <typename Types, size_t... ISeq> class par_iterator {
public:
  template <size_t I>
  using slice = iterator<typename Types::template ith_type<I>>;

  explicit par_iterator(slice<ISeq> const &...slices) : slices{slices...} {};

  at_expr<Types> operator*() const { return {*slices.template get<ISeq>()...}; }

  par_iterator &operator++() {
    (..., ++slices.template get<ISeq>());
    return *this;
  }

  par_iterator operator++(int) {
    return par_iterator(slices.template get<ISeq>()++...);
  }

  template <typename Idx> par_iterator operator+(Idx const &offset) const {
    return par_iterator(slices.template get<ISeq>() + offset...);
  }

  template <typename Idx> par_iterator &operator+=(Idx const &offset) {
    (..., advance<ISeq>(offset));
    return *this;
  }

  par_iterator &operator--() {
    (..., --slices.template get<ISeq>());
    return *this;
  }

  par_iterator operator--(int) {
    return par_iterator(slices.template get<ISeq>()--...);
  }

  template <typename Idx> par_iterator operator-(Idx const &offset) const {
    return par_iterator(slices.template get<ISeq>() - offset...);
  }

  std::ptrdiff_t operator-(par_iterator<Types, ISeq...> const &other) {
    return slices.template get<0>() - other.slices.template get<0>();
  }

  template <typename Idx> par_iterator &operator-=(Idx const &offset) {
    (..., backtrack<ISeq>(offset));
    return *this;
  }

  bool operator<(par_iterator const &other) const {
    return slices.template get<0>() < other.slices.template get<0>();
  }

  bool operator<=(par_iterator const &other) const {
    return slices.template get<0>() <= other.slices.template get<0>();
  }

  bool operator>(par_iterator const &other) const {
    return slices.template get<0>() > other.slices.template get<0>();
  }

  bool operator>=(par_iterator const &other) const {
    return slices.template get<0>() >= other.slices.template get<0>();
  }

  bool operator==(par_iterator const &other) const {
    return slices.template get<0>() == other.slices.template get<0>();
  }

  bool operator!=(par_iterator const &other) const {
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

namespace std {
template <typename Types, size_t... ISeq>
struct iterator_traits<nitro::par_iterator<Types, ISeq...>> {
  using value_type = Types;
  using difference_type = std::ptrdiff_t;
  using reference = nitro::par_at_expr<Types, ISeq...>;
  using iterator_category = std::random_access_iterator_tag;
};
} // namespace std