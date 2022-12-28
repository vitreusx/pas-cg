#pragma once
#include "../cuda_interop.h"
#include "const_iterator.h"
#include "lane.h"

namespace nitro::def {
template <typename T>
class const_view {
public:
  const_view() = default;
  explicit const_view(T const *data, int n) : data{data}, n{n} {};

  __host__ __device__ T const &operator[](int idx) const {
    return data[idx];
  }

  template <typename Idxes, typename = std::enable_if_t<is_lane_like_v<Idxes>>>
  auto operator[](Idxes idxes) const {
    using Data = lane<T, lane_size_v<Idxes>, lane_width_v<Idxes>>;
    return gather<Data>(data, idxes);
  }

  template <typename Idxes, typename Mask,
            typename =
                std::enable_if_t<is_lane_like_v<Idxes> && is_lane_like_v<Mask>>>
  auto operator[](std::pair<Idxes, Mask> idxes_mask) const {
    using Data = lane<T, lane_size_v<Idxes>, lane_width_v<Idxes>>;
    auto const &[idxes, mask] = idxes_mask;
    return masked_gather<Data>(data, idxes, mask);
  }

  template <typename Idx>
  __host__ __device__ decltype(auto) at(Idx idx) const {
    return (*this)[idx];
  }

  template <std::size_t N, std::size_t W = opt_width_v>
  lane<T, N, W> at_lane(int idx) const {
    return construct<lane<T, N, W>>(data + N * idx);
  }

  __host__ __device__ int size() const {
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