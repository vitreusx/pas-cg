#pragma once
#include "const_iterator.h"

namespace nitro::debug {
template <typename T> class const_view {
public:
  const_view() = default;
  explicit const_view(std::deque<T> const *data) : data{data} {
    if (data)
      n = data->size();
    else
      n = 0;
  };

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

  decltype(auto) front() const {
    return data->front();
  }

  decltype(auto) back() const {
    return data->back();
  }

protected:
  std::deque<T> const *data = nullptr;
  int n = 0;
};
} // namespace nitro::debug