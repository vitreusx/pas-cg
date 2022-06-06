#pragma once
#include <iterator>

namespace nitro::aos {
template <typename T> class iterator {
public:
  explicit iterator(T *ptr) : ptr{ptr} {};

  T &operator*() const {
    return *ptr;
  }

  auto &operator++() {
    ++ptr;
    return *this;
  }

  auto operator++(int) {
    return iterator(ptr++);
  }

  template <typename Idx> auto operator+(Idx const &offset) const {
    return iterator(ptr + offset);
  }

  template <typename Idx> auto &operator+=(Idx const &offset) {
    ptr += offset;
    return *this;
  }

  auto &operator--() {
    --ptr;
    return *this;
  }

  auto operator--(int) {
    return iterator(ptr--);
  }

  template <typename Idx> auto operator-(Idx const &offset) const {
    return iterator(ptr - offset);
  }

  template <typename Idx> auto &operator-=(Idx const &offset) {
    ptr -= offset;
    return *this;
  }

  std::ptrdiff_t operator-(iterator<T> const &other) {
    return ptr - other.ptr;
  }

  bool operator<(iterator const &other) const {
    return ptr < other.ptr;
  }

  bool operator<=(iterator const &other) const {
    return ptr <= other.ptr;
  }

  bool operator>(iterator const &other) const {
    return ptr > other.ptr;
  }

  bool operator>=(iterator const &other) const {
    return ptr >= other.ptr;
  }

  bool operator==(iterator const &other) const {
    return ptr == other.ptr;
  }

  bool operator!=(iterator const &other) const {
    return ptr != other.ptr;
  }

private:
  T *ptr;
};
} // namespace nitro::aos

namespace std {
template <typename T> struct iterator_traits<nitro::aos::iterator<T>> {
  using value_type = T;
  using difference_type = std::ptrdiff_t;
  using reference = T &;
  using iterator_category = std::random_access_iterator_tag;
};
} // namespace std