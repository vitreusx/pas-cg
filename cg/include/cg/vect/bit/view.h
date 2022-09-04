#pragma once
#include "const_view.h"
#include "iterator.h"
#include "lane_ref.h"

namespace nitro::bit {
class view {
public:
  inline view() = default;
  inline explicit view(byte *data, int n) : data{data}, n{n} {};

  inline ref operator[](int idx) const {
    auto mask = (byte)1 << (idx % num_bits);
    return ref(data + idx / num_bits, mask);
  }

  inline ref at(int idx) const {
    return ref(data + byte_offset(idx), byte_mask(idx));
  }

  template <std::size_t N, std::size_t W>
  inline lane_ref<N, W> at_lane(int idx) const {
    return lane_ref<N, W>(data + byte_offset(N * idx));
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

  inline operator const_view() const {
    return const_view(data, n);
  }

  inline auto begin() const {
    return iterator(data, 0);
  }

  inline auto end() const {
    return iterator(data, n);
  }

private:
  byte *data = nullptr;
  int n = 0;
};
} // namespace nitro::bit