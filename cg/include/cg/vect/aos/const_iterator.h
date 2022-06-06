#pragma once
#include <memory>

namespace nitro::aos {
template <typename T> class const_iterator {
public:
  explicit const_iterator(T const *ptr) : ptr{ptr} {};

  T const &operator*() const {
    return *ptr;
  }

  auto &operator++() {
    ++ptr;
    return *this;
  }

  auto operator++(int) {
    return const_iterator(ptr++);
  }

  template <typename Idx> auto operator+(Idx const &offset) const {
    return def_const_iterator(ptr + offset);
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
    return const_iterator(ptr--);
  }

  template <typename Idx> auto operator-(Idx const &offset) const {
    return const_iterator(ptr - offset);
  }

  template <typename Idx> auto &operator-=(Idx const &offset) {
    ptr -= offset;
    return *this;
  }

  std::ptrdiff_t operator-(const_iterator<T> const &other) {
    return ptr - other.ptr;
  }

  bool operator<(const_iterator const &other) const {
    return ptr < other.ptr;
  }

  bool operator<=(const_iterator const &other) const {
    return ptr <= other.ptr;
  }

  bool operator>(const_iterator const &other) const {
    return ptr > other.ptr;
  }

  bool operator>=(const_iterator const &other) const {
    return ptr >= other.ptr;
  }

  bool operator==(const_iterator const &other) const {
    return ptr == other.ptr;
  }

  bool operator!=(const_iterator const &other) const {
    return ptr != other.ptr;
  }

private:
  T const *ptr;
};
} // namespace nitro::aos

namespace std {
template <typename T> struct iterator_traits<nitro::aos::const_iterator<T>> {
  using value_type = T const;
  using difference_type = std::ptrdiff_t;
  using reference = T const &;
  using iterator_category = std::random_access_iterator_tag;
};
} // namespace std