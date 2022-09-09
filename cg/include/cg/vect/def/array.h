#pragma once
#include "const_iterator.h"
#include "const_view.h"
#include "iterator.h"
#include "view.h"

namespace nitro::def {
template <typename T, std::size_t N>
class array {
public:
  template <typename... Args>
  array(Args &&...args) : data{std::forward<Args>(args)...} {}

  int size() const {
    return N;
  }

  decltype(auto) operator[](int idx) {
    return data[idx];
  }

  template <typename Idxes, typename = std::enable_if_t<is_lane_like_v<Idxes>>>
  auto operator[](Idxes idxes) {
    return sparse_ref<T, Idxes>(data.data(), idxes);
  }

  template <typename Idxes, typename Mask,
            typename =
                std::enable_if_t<is_lane_like_v<Idxes> && is_lane_like_v<Mask>>>
  auto operator[](std::pair<Idxes, Mask> idxes_mask) {
    auto const &[idxes, mask] = idxes_mask;
    return masked_ref(data.data(), idxes, mask);
  }

  decltype(auto) operator[](int idx) const {
    return data[idx];
  }

  template <typename Idxes, typename = std::enable_if_t<is_lane_like_v<Idxes>>>
  auto operator[](Idxes idxes) const {
    using Data = lane<T, lane_size_v<Idxes>, lane_width_v<Idxes>>;
    return gather<Data>(data.data(), idxes);
  }

  template <typename Idxes, typename Mask,
            typename =
                std::enable_if_t<is_lane_like_v<Idxes> && is_lane_like_v<Mask>>>
  auto operator[](std::pair<Idxes, Mask> idxes_mask) const {
    using Data = lane<T, lane_size_v<Idxes>, lane_width_v<Idxes>>;
    auto const &[idxes, mask] = idxes_mask;
    return masked_gather<Data>(data.data(), idxes, mask);
  }

  template <typename Idx>
  decltype(auto) at(Idx idx) {
    return (*this)[idx];
  }

  template <typename Idx>
  decltype(auto) at(Idx idx) const {
    return (*this)[idx];
  }

  template <std::size_t N_, std::size_t W = opt_width_v>
  lane_ref<T, N_, W> at_lane(int idx) {
    return lane_ref<T, N_, W>(data.data() + N_ * idx);
  }

  template <std::size_t N_, std::size_t W = opt_width_v>
  lane<T, N_, W> at_lane(int idx) const {
    return construct<lane<T, N_, W>>(data.data() + N_ * idx);
  }

  template <std::size_t N_>
  int num_lanes() const {
    return size() / N_;
  }

  template <std::size_t N_>
  int final_idx() const {
    return N_ * num_lanes<N_>();
  }

  void clear() {
    std::destroy_n(data.data(), size());
  }

  auto begin() {
    return iterator<T>(data.data());
  }

  auto begin() const {
    return const_iterator<T>(data.data());
  }

  auto end() {
    return iterator<T>(data.data() + size());
  }

  auto end() const {
    return const_iterator<T>(data.data() + size());
  }

  operator view<T>() {
    return view<T>(data.data(), size());
  }

  operator const_view<T>() const {
    return const_view<T>(data.data(), size());
  }

private:
  std::array<T, N> data;
};
} // namespace nitro::def