#pragma once
#include "../expr.h"
#include "../tuple.h"
#include "decl.h"
#include "def.h"
#include <cstddef>
#include <utility>

namespace nitro {

template <typename Types, size_t... ISeq>
class par_allocator
    : public tuple<allocator<typename Types::template ith_type<ISeq>>...> {
public:
  template <size_t I>
  using slice = allocator<typename Types::template ith_type<I>>;

  using Base = tuple<slice<ISeq>...>;
  using Base::Base;

  par_allocator() : Base{slice<ISeq>()...} {};
};
} // namespace nitro