#pragma once
#include "lane.h"

namespace nitro::def {
template <typename T, std::size_t N, std::size_t W>
class lane_ref {
public:
  using Data = lane<T, N, W>;

  explicit lane_ref(T *p) : p{p}, data{construct<Data>(p)} {}

  operator Data const &() const {
    return data;
  }

  template <typename U>
  auto &operator=(U const &value) {
    data = value;
    store(data, p);
    return *this;
  }

  template <typename U>
  friend void store(lane_ref const &ref, U *dst) {
    store(static_cast<Data const &>(ref), dst);
  }

private:
  T *p;
  Data data;
};

template <typename T, std::size_t N, std::size_t W>
struct is_lane_like<lane_ref<T, N, W>> : std::true_type {};

template <typename T, std::size_t N, std::size_t W>
struct lane_type<lane_ref<T, N, W>> {
  using type = lane_type_t<typename lane_ref<T, N, W>::Data>;
};

template <typename T, std::size_t N, std::size_t W>
struct lane_size<lane_ref<T, N, W>> {
  static constexpr std::size_t value =
      lane_size_v<typename lane_ref<T, N, W>::Data>;
};

} // namespace nitro::def