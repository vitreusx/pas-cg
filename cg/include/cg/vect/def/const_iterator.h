#pragma once
#include <memory>

namespace nitro::def {
template <typename T>
class const_iterator {
public:
  explicit const_iterator(T const *ptr) : ptr{ptr} {};
  explicit const_iterator(T const &ref) : ptr{&ref} {}

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

  auto operator+(int offset) const {
    return def_const_iterator(ptr + offset);
  }

  auto &operator+=(int offset) {
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

  auto operator-(int offset) const {
    return const_iterator(ptr - offset);
  }

  auto &operator-=(int offset) {
    ptr -= offset;
    return *this;
  }

  std::ptrdiff_t operator-(const_iterator<T> const &other) const {
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
} // namespace nitro::def

namespace std {
template <typename T>
struct iterator_traits<nitro::def::const_iterator<T>> {
  using value_type = T const;
  using difference_type = std::ptrdiff_t;
  using reference = T const &;
  using iterator_category = std::random_access_iterator_tag;
};
} // namespace std