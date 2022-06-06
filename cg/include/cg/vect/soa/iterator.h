#pragma once
#include "../aos/iterator.h"
#include "../common/is_indexed.h"
#include "../common/tuple.h"

namespace nitro::soa {
template <bool IsIndexed, typename T>
struct iterator_impl;

template <typename T>
using iterator = typename iterator_impl<is_indexed_v<T>, T>::type;

template <typename T>
struct iterator_impl<false, T> {
  using type = aos::iterator<T>;
};

template<typename T>
class soa_iterator;

template <typename... Types> class soa_iterator<type_list<Types...>> {
public:
  explicit soa_iterator(iterator<Types> const &...slices) : slices{slices...} {};

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

  tuple<iterator<Types>...> slices;
};
} // namespace nitro::soa