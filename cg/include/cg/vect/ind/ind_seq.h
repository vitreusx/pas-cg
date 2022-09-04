#pragma once
#include <utility>

namespace nitro::ind {
template <std::size_t... Idxes>
using ind_seq = std::index_sequence<Idxes...>;
}