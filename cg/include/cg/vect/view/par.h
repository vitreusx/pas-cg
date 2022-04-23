#pragma once
#include "../tuple.h"
#include "decl.h"
#include "../iterator.h"
#include "../const_iterator.h"
#include <cstddef>

namespace nitro {
template <typename Types, typename Idx, size_t... ISeq> class par_view {
public:
  template <size_t I>
  using slice = view<typename Types::template ith_type<I>, Idx>;

  par_view() : slices{slice<ISeq>()...} {};
  explicit par_view(slice<ISeq> const &...slices) : slices{slices...} {};

  at_expr<Types> operator[](Idx idx) const {
    return {slices.template get<ISeq>()[idx]...};
  }

  at_expr<Types> at(Idx idx) const { return (*this)[idx]; }

  template <size_t N> lane_at_expr<Types, N> lane_at(Idx idx) const {
    return {slices.template get<ISeq>().template lane_at<N>(idx)...};
  }

  template <size_t N> Idx num_lanes() const { return size() / N; }

  Idx size() const { return slices.template get<0>().size(); }

  operator par_const_view<Types, Idx, ISeq...>() const {
    return par_const_view<Types, Idx, ISeq...>(slices.template get<ISeq>()...);
  }

  iterator<Types> begin() {
    return par_iterator<Types, ISeq...>(slices.template get<ISeq>().begin()...);
  }

  const_iterator<Types> begin() const {
    return par_const_iterator<Types, ISeq...>(slices.template get<ISeq>().begin()...);
  }

  iterator<Types> end() {
    return par_iterator<Types, ISeq...>(slices.template get<ISeq>().end()...);
  }

  const_iterator<Types> end() const {
    return par_const_iterator<Types, ISeq...>(slices.template get<ISeq>().end()...);
  }

private:
  tuple<slice<ISeq>...> slices;
};
} // namespace nitro