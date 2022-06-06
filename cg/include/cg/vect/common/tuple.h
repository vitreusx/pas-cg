#pragma once
#include "type_list.h"
#include <utility>

namespace nitro {
template <typename... Types>
class tuple;

template <>
class tuple<> {
public:
  static auto Idxes() {
    return std::make_index_sequence<0>{};
  }

  static auto Types() {
    return type_list<>{};
  }
};

template <typename Head, typename... Tail>
class tuple<Head, Tail...> {
public:
  template <typename HeadArg, typename... TailArgs>
  tuple(HeadArg &&head, TailArgs &&...tail)
      : head{std::forward<HeadArg>(head)}, tail{std::forward<TailArgs>(
                                               tail)...} {};

  template <typename E>
  static tuple<Head, Tail...> from_expr(E const &e) {
    return _from_expr(e, e.Idxes());
  }

  template <std::size_t I>
  decltype(auto) get() {
    if constexpr (I == 0)
      return (head);
    else
      return tail.template get<I - 1>();
  }

  template <std::size_t I>
  decltype(auto) get() const {
    if constexpr (I == 0)
      return (head);
    else
      return tail.template get<I - 1>();
  }

  static auto Idxes() {
    return std::make_index_sequence<1 + sizeof...(Tail)>{};
  }

  static auto Types() {
    return type_list<Head, Tail...>{};
  }

private:
  Head head;
  tuple<Tail...> tail;

  template <typename E, std::size_t... Idxes>
  static tuple<Head, Tail...> _from_expr(E const &e,
                                         std::index_sequence<Idxes...>) {
    return {e.template get<Idxes>()...};
  }
};
} // namespace nitro