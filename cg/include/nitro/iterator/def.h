#pragma once
#include "../at.h"
#include "decl.h"
#include <iterator>

namespace nitro {
template <typename T> class def_iterator {
public:
  explicit def_iterator(T *ptr) : ptr{ptr} {};

  at_expr<T> operator*() const { return *ptr; }

  def_iterator &operator++() {
    ++ptr;
    return *this;
  }

  def_iterator operator++(int) { return def_iterator(ptr++); }

  template <typename Idx> def_iterator operator+(Idx const &offset) const {
    return def_iterator(ptr + offset);
  }

  template <typename Idx> def_iterator &operator+=(Idx const &offset) {
    ptr += offset;
    return *this;
  }

  def_iterator &operator--() {
    --ptr;
    return *this;
  }

  def_iterator operator--(int) { return def_iterator(ptr--); }

  template <typename Idx> def_iterator operator-(Idx const &offset) const {
    return def_iterator(ptr - offset);
  }

  template <typename Idx> def_iterator &operator-=(Idx const &offset) {
    ptr -= offset;
    return *this;
  }

  std::ptrdiff_t operator-(def_iterator<T> const &other) {
    return ptr - other.ptr;
  }

  bool operator<(def_iterator const &other) const { return ptr < other.ptr; }

  bool operator<=(def_iterator const &other) const { return ptr <= other.ptr; }

  bool operator>(def_iterator const &other) const { return ptr > other.ptr; }

  bool operator>=(def_iterator const &other) const { return ptr >= other.ptr; }

  bool operator==(def_iterator const &other) const { return ptr == other.ptr; }

  bool operator!=(def_iterator const &other) const { return ptr != other.ptr; }

private:
  T *ptr;
};
} // namespace nitro

namespace std {
template <typename T> struct iterator_traits<nitro::def_iterator<T>> {
  using value_type = T;
  using difference_type = std::ptrdiff_t;
  using reference = nitro::at_expr<T>;
  using iterator_category = std::random_access_iterator_tag;
};
} // namespace std