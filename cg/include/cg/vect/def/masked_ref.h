#pragma once
#include "lane.h"

namespace nitro::def {
template <typename T, typename Idxes, typename Mask>
class masked_ref {
public:
  using Data = lane<T, lane_size_v<Idxes>, lane_width_v<Idxes>>;

  explicit masked_ref(T *p, Idxes const &idxes, Mask const &mask)
      : p{p}, idxes{idxes}, mask{mask} {}

  operator Data() const {
    return masked_gather<Data>(p, idxes, mask);
  }

  template <typename U>
  auto &operator=(U const &value) {
    auto data = masked_gather<Data>(p, idxes, mask);
    data = value;
    masked_scatter(data, p, idxes, mask);
    return *this;
  }

  template <typename U>
  friend void store(masked_ref const &ref, U *dst) {
    store(static_cast<Data const &>(ref), dst);
  }

private:
  T *p;
  Idxes idxes;
  Mask mask;
};

template <typename T, typename Idxes, typename Mask>
struct is_lane_like<masked_ref<T, Idxes, Mask>> : std::true_type {};

template <typename T, typename Idxes, typename Mask>
struct lane_type<masked_ref<T, Idxes, Mask>> {
  using type = lane_type_t<typename masked_ref<T, Idxes, Mask>::Data>;
};

template <typename T, typename Idxes, typename Mask>
struct lane_size<masked_ref<T, Idxes, Mask>> {
  static constexpr std::size_t value =
      lane_size_v<typename masked_ref<T, Idxes, Mask>::Data>;
};
} // namespace nitro::def