#pragma once
#include "lane.h"

namespace nitro::def {
template <typename T, typename Idxes>
class sparse_ref {
public:
  using Data = lane<T, lane_size_v<Idxes>, lane_width_v<Idxes>>;

  explicit sparse_ref(T *p, Idxes const &idxes) : p{p}, idxes{idxes} {}

  operator Data() const {
    return gather<Data>(p, idxes);
  }

  template <typename U>
  auto &operator=(U const &value) {
    auto data = gather<Data>(p, idxes);
    data = value;
    scatter(data, p, idxes);
    return *this;
  }

  template <typename U>
  friend void store(sparse_ref const &ref, U *dst) {
    store(static_cast<Data const &>(ref), dst);
  }

private:
  T *p;
  Idxes idxes;
};

template <typename T, typename Idxes>
struct is_lane_like<sparse_ref<T, Idxes>> : std::true_type {};

template <typename T, typename Idxes>
struct lane_type<sparse_ref<T, Idxes>> {
  using type = lane_type_t<typename sparse_ref<T, Idxes>::Data>;
};

template <typename T, typename Idxes>
struct lane_size<sparse_ref<T, Idxes>> {
  static constexpr std::size_t value =
      lane_size_v<typename sparse_ref<T, Idxes>::Data>;
};
} // namespace nitro::def