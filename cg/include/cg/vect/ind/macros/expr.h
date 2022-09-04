#pragma once
#include "../ind_seq.h"
#include "expr_get_body.h"
#include "iexpr_named.h"
#include "vfunc.h"

#define EXPR(...)                                                              \
  _IEXPR_NAMED(__VA_ARGS__)                                                    \
                                                                               \
  template <std::size_t I>                                                     \
  decltype(auto) get() {                                                       \
    _EXPR_GET_BODY(, __VA_ARGS__)                                              \
  }                                                                            \
                                                                               \
  template <std::size_t I>                                                     \
  decltype(auto) get() const {                                                 \
    _EXPR_GET_BODY(const, __VA_ARGS__)                                         \
  }                                                                            \
                                                                               \
  static auto Idxes() {                                                        \
    return std::make_index_sequence<NARG(__VA_ARGS__)>{};                      \
  }
