#pragma once
#include "const_iterator.h"

namespace nitro::aos {
template <typename T> class const_view {
public:
  const_view() = default;
  explicit const_view(T const *data, int n) : data{data}, n{n} {};

  T const &operator[](int idx) const {
    return data[idx];
  }

  T const &at(int idx) const {
    return (*this)[idx];
  }

  int size() const {
    return n;
  }

  auto begin() const {
    return const_iterator(data);
  }

  auto end() const {
    return const_iterator(data + n);
  }

protected:
  T const *data = nullptr;
  int n = 0;
};
} // namespace nitro::aos