#pragma once
#include "../cuda_interop.h"
#include "const_view.h"
#include "iterator.h"
#include "lane_ref.h"
#include "masked_ref.h"
#include "sparse_ref.h"

namespace nitro::def {
template <typename T>
class view {
public:
  view() = default;
  explicit view(T *data, int n) : data{data}, n{n} {};

  __host__ __device__ T &operator[](int idx) const {
    return data[idx];
  }

  template <typename Idxes, typename = std::enable_if_t<is_lane_like_v<Idxes>>>
  auto operator[](Idxes idxes) const {
    return sparse_ref(data, idxes);
  }

  template <typename Idxes, typename Mask,
            typename =
                std::enable_if_t<is_lane_like_v<Idxes> && is_lane_like_v<Mask>>>
  auto operator[](std::pair<Idxes, Mask> idxes_mask) const {
    auto const &[idxes, mask] = idxes_mask;
    return masked_ref(data, idxes, mask);
  }

  template <typename Idx>
  __host__ __device__ decltype(auto) at(Idx idx) const {
    return (*this)[idx];
  }

  template <std::size_t N, std::size_t W = opt_width_v>
  lane_ref<T, N, W> at_lane(int idx) const {
    return lane_ref<T, N, W>(data + N * idx);
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