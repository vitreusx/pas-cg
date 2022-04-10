#pragma once
#include <cg/nl/exclusion.h>
#include <cg/nl/pair.h>

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

template <typename E1, typename E2>
bool operator<(pair_expr<E1> const& pair1, pair_expr<E2> const& pair2) {
  return std::make_pair(pair1.i1(), pair1.i2()) < std::make_pair(pair2.i1(), pair2.i2());
}

} // namespace cg::nl