#pragma once
#include "config.h"
#include "type_list.h"
#include <cstddef>
#include <utility>

namespace nitro {
template <typename E> struct ind_expr {
  template <size_t I>  decltype(auto) get() {
    return static_cast<E &>(*this).template get<I>();
  }

  template <size_t I>  decltype(auto) get() const {
    return static_cast<E const &>(*this).template get<I>();
  }
};

template <typename T> struct is_indexed_impl : std::false_type {};

template <typename T> constexpr bool is_indexed = is_indexed_impl<T>::value;

} // namespace nitro