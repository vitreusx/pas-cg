#pragma once
#include "lane.h"

namespace nitro::def {
template <typename T, typename Idxes>
using sparse_const_ref = lane<T, lane_size_v<Idxes>, lane_width_v<Idxes>>;

template <typename T, typename Idxes, typename Mask>
using masked_const_ref = lane<T, lane_size_v<Idxes>, lane_width_v<Idxes>>;
} // namespace nitro::def