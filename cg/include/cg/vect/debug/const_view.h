#pragma once
#include "const_iterator.h"

namespace nitro::debug {
template <typename T> class const_view {
public:
  const_view() = default;
  explicit const_view(std::deque<T> const *data) : data{data} {};

  T const &operator[](int idx) const {
    return data->at(idx);
  }

  T const &at(int idx) const {
    return (*this)[idx];
  }

  int size() const {
    if (!data)
      return 0;
    else
      return data->size();
  }

  auto begin() const {
    return data->begin();
  }

  auto end() const {
    return data->end();
  }

protected:
  std::deque<T> const *data = nullptr;
};
} // namespace nitro::debug