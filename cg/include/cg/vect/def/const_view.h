#pragma once
#include "const_iterator.h"
#include "lane_const_ref.h"

namespace nitro::def {
template <typename T>
class const_view {
public:
  const_view() = default;
  explicit const_view(T const *data, int n) : data{data}, n{n} {};

  T const &operator[](int idx) const {
    return data[idx];
  }

  T const &at(int idx) const {
    return (*this)[idx];
  }

  template <std::size_t N, std::size_t W>
  lane_const_ref<T, N, W> at_lane(int idx) const {
    return lane_const_ref<T, N, W>(data + N * idx);
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

  auto begin() const {
    return const_iterator<T>(data);
  }

  auto end() const {
    return const_iterator<T>(data + n);
  }

protected:
  T const *data = nullptr;
  int n = 0;
};
} // namespace nitro::def