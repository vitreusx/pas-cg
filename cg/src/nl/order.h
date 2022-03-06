#pragma once
#include "nl/exclusion.h"
#include "nl/pair.h"

namespace cg::nl {
template <typename E1, typename E2>
bool operator<(exclusion_expr<E1> const &excl, pair_expr<E2> const &pair) {
  return std::make_pair(excl.i1(), excl.i1()) <
         std::make_pair(pair.i1(), pair.i2());
}

template <typename E1, typename E2>
bool operator>(exclusion_expr<E1> const &excl, pair_expr<E2> const &pair) {
  return std::make_pair(excl.i1(), excl.i1()) >
         std::make_pair(pair.i1(), pair.i2());
}
} // namespace cg::nl