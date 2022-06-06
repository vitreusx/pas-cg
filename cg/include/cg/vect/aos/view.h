#pragma once
#include "const_view.h"
#include "iterator.h"

namespace nitro::aos {
template <typename T> class view {
public:
  view() = default;
  explicit view(T *data, int n) : data{data}, n{n} {};

  T &operator[](int idx) const {
    return data[idx];
  }

  T &at(int idx) const {
    return (*this)[idx];
  }

  int size() const {
    return n;
  }

  operator const_view<T>() const {
    return const_view<T>(data, n);
  }

  auto begin() {
    return iterator(data);
  }

  auto begin() const {
    return const_iterator(data);
  }

  auto end() {
    return iterator(data + n);
  }

  auto end() const {
    return const_iterator(data + n);
  }

private:
  T *data = nullptr;
  int n = 0;
};
} // namespace nitro::aos