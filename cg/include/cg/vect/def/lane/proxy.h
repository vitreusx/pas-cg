#pragma once
#include "../repr.h"
#include "vcl.h"

namespace nitro::def {
template <typename T, std::size_t N, std::size_t W>
using proxy_lane = vcl_lane<repr_t<T>, N, W>;
} // namespace nitro::def