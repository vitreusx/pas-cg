#pragma once
#include <iterator>

namespace nitro::def {
template <typename T>
class iterator {
public:
  explicit iterator(T *ptr) : ptr{ptr} {}
  explicit iterator(T &ref) : ptr{&ref} {}

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

  auto operator+(int offset) const {
    return iterator(ptr + offset);
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
    return iterator(ptr--);
  }

  auto operator-(int offset) const {
    return iterator(ptr - offset);
  }

  auto &operator-=(int offset) {
    ptr -= offset;
    return *this;
  }

  std::ptrdiff_t operator-(iterator<T> const &other) const {
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
} // namespace nitro::def

namespace std {
template <typename T>
struct iterator_traits<nitro::def::iterator<T>> {
  using value_type = T;
  using difference_type = std::ptrdiff_t;
  using reference = T &;
  using iterator_category = std::random_access_iterator_tag;
};
} // namespace std