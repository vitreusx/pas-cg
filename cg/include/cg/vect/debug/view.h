#pragma once
#include "const_view.h"
#include "iterator.h"

namespace nitro::debug {
template <typename T> class view {
public:
  view() = default;
  explicit view(std::deque<T> *data) : data{data} {
    if (data)
      n = data->size();
    else
      n = 0;
  };

  T &operator[](int idx) const {
    return data->at(idx);
  }

  T &at(int idx) const {
    return (*this)[idx];
  }

  int size() const {
    if (!data)
      return 0;
    else
      return data->size();
  }

  operator const_view<T>() const {
    return const_view<T>(data);
  }

  auto begin() {
    return data->begin();
  }

  auto begin() const {
    return data->cbegin();
  }

  auto end() {
    return data->end();
  }

  auto end() const {
    return data->cend();
  }

  decltype(auto) front() const {
    return data->front();
  }

  decltype(auto) back() const {
    return data->back();
  }

private:
  std::deque<T> *data = nullptr;
  int n = 0;
};
} // namespace nitro::debug