#pragma once
#include "const_iterator.h"
#include "lane_const_ref.h"

namespace nitro::bit {
class const_view {
public:
  inline const_view() = default;
  inline explicit const_view(byte const *data, int n) : data{data}, n{n} {};

  inline const_ref operator[](int idx) const {
    auto mask = (byte)1 << (idx % num_bits);
    return const_ref(data + idx / num_bits, mask);
  }

  inline const_ref at(int idx) const {
    return const_ref(data + byte_offset(idx), byte_mask(idx));
  }

  template <std::size_t N, std::size_t W>
  lane_const_ref<N, W> at_lane(int idx) const {
    return lane_const_ref<N, W>(data + byte_offset(N * idx));
  }

  inline int size() const {
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

  inline auto begin() const {
    return const_iterator(data);
  }

  inline auto end() const {
    return const_iterator(data + n);
  }

protected:
  byte const *data = nullptr;
  int n = 0;
};
} // namespace nitro::bit