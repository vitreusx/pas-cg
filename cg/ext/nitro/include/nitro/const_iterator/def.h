#pragma once
#include "../const_at.h"
#include "decl.h"
#include <iterator>

namespace nitro {
template <typename T> class def_const_iterator {
public:
  explicit def_const_iterator(T const *ptr) : ptr{ptr} {};

  const_at_expr<T> operator*() const { return *ptr; }

  def_const_iterator &operator++() {
    ++ptr;
    return *this;
  }

  def_const_iterator operator++(int) { return def_const_iterator(ptr++); }

  template <typename Idx>
  def_const_iterator operator+(Idx const &offset) const {
    return def_const_iterator(ptr + offset);
  }

  template <typename Idx> def_const_iterator &operator+=(Idx const &offset) {
    ptr += offset;
    return *this;
  }

  def_const_iterator &operator--() {
    --ptr;
    return *this;
  }

  def_const_iterator operator--(int) { return def_const_iterator(ptr--); }

  template <typename Idx>
  def_const_iterator operator-(Idx const &offset) const {
    return def_const_iterator(ptr - offset);
  }

  template <typename Idx> def_const_iterator &operator-=(Idx const &offset) {
    ptr -= offset;
    return *this;
  }

  std::ptrdiff_t operator-(def_const_iterator<T> const &other) {
    return ptr - other.ptr;
  }

  bool operator<(def_const_iterator const &other) const {
    return ptr < other.ptr;
  }

  bool operator<=(def_const_iterator const &other) const {
    return ptr <= other.ptr;
  }

  bool operator>(def_const_iterator const &other) const {
    return ptr > other.ptr;
  }

  bool operator>=(def_const_iterator const &other) const {
    return ptr >= other.ptr;
  }

  bool operator==(def_const_iterator const &other) const {
    return ptr == other.ptr;
  }

  bool operator!=(def_const_iterator const &other) const {
    return ptr != other.ptr;
  }

private:
  T const* ptr;
};
} // namespace nitro

namespace std {
template <typename T> struct iterator_traits<nitro::def_const_iterator<T>> {
  using value_type = T const;
  using difference_type = std::ptrdiff_t;
  using reference = nitro::const_at_expr<T>;
  using iterator_category = std::random_access_iterator_tag;
};
} // namespace std