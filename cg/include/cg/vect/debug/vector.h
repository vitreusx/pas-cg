#pragma once
#include "const_view.h"
#include "view.h"
#include <vector>

namespace nitro::debug {
template <typename T, typename Alloc = std::allocator<T>> class vector {
public:
  vector() = default;
  explicit vector(Alloc alloc) : base{alloc} {}
  explicit vector(int n, T const &init = T(), Alloc alloc = Alloc())
      : base(n, init, alloc) {}

  int size() const {
    return base.size();
  }

  int capacity() const {
    return base.capacity();
  }

  T &operator[](int idx) {
    return at(idx);
  }

  T &at(int idx) {
    return base.at(idx);
  }

  T const &operator[](int idx) const {
    return at(idx);
  }

  T const &at(int idx) const {
    return base.at(idx);
  }

  void clear() {
    base.clear();
  }

  void reserve(int new_capacity) {
    base.reserve(new_capacity);
  }

  void resize(int new_size, T const &init = T()) {
    base.resize(new_size, init);
  }

  void shrink(int new_size) {
    base.erase(base.begin() + new_size, base.end());
  }

  void push_back(T const &value) {
    base.push_back(value);
  }

  template <typename... Args> T &emplace_back(Args &&...args) {
    return base.emplace_back(std::forward<Args>(args)...);
  }

  auto begin() {
    return base.begin();
  }

  auto end() {
    return base.end();
  }

  auto begin() const {
    return base.begin();
  }

  auto end() const {
    return base.end();
  }

  operator view<T>() {
    return view<T>(&base);
  }

  operator const_view<T>() const {
    return const_view<T>(&base);
  }

private:
  std::deque<T, Alloc> base;
};
} // namespace nitro::debug