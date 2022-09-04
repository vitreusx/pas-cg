#pragma once
#include "const_view.h"
#include "iterator.h"
#include "lane_ref.h"

namespace nitro::def {
template <typename T>
class view {
public:
  view() = default;
  explicit view(T *data, int n) : data{data}, n{n} {};

  T &operator[](int idx) const {
    return data[idx];
  }

  T &at(int idx) const {
    return (*this)[idx];
  }

  template <std::size_t N, std::size_t W>
  lane_ref<T, N, W> at_lane(int idx) const {
    return lane_ref<T, N, W>(data + N * idx);
  }

  int size() const {
    return n;
  }

  template <std::size_t N>
  int num_lanes() const {
    return size() / N;
  }

  template <std::size_t N>
  int final_idx() const {
    return N * num_lanes<N>();
  }

  operator const_view<T>() const {
    return const_view<T>(data, n);
  }

  auto begin() const {
    return iterator<T>(data);
  }

  auto end() const {
    return iterator<T>(data + n);
  }

private:
  T *data = nullptr;
  int n = 0;
};
} // namespace nitro::def