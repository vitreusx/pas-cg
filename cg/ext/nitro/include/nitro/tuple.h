#pragma once
#include "expr.h"
#include "indexed.h"
#include "type_list.h"
#include <utility>

namespace nitro {
template <typename... Types> class tuple;

template <> class tuple<> : public ind_expr<tuple<>> {
public:
  template <size_t I> using ith_type = void;
  static constexpr size_t num_types = 0;
};

template <typename Head, typename... Tail>
class tuple<Head, Tail...> : public ind_expr<tuple<Head, Tail...>> {
public:
   tuple(Head const &head, Tail const &...tail)
      : head{head}, tail{tail...} {}

  template <typename E>
   tuple(ind_expr<E> const &e)
      : tuple(e, std::make_index_sequence<1 + sizeof...(Tail)>{}) {}

  template <size_t I>  decltype(auto) get() {
    if constexpr (I == 0)
      return (head);
    else
      return tail.template get<I - 1>();
  }

  template <size_t I>  decltype(auto) get() const {
    if constexpr (I == 0)
      return (head);
    else
      return tail.template get<I - 1>();
  }

  static constexpr size_t num_types = type_list<Head, Tail...>::num_types;

  template <size_t I>
  using ith_type = typename type_list<Head, Tail...>::template ith_type<I>;

private:
  Head head;
  tuple<Tail...> tail;

  template <typename E, size_t... ISeq>
   tuple(ind_expr<E> const &e, std::index_sequence<ISeq...>)
      : tuple{e.template get<ISeq>()...} {}
};

template <typename E, typename... Types> struct expr_impl<E, tuple<Types...>> {
  using type = ind_expr<E>;
};

template <typename E, typename... Types>
struct auto_expr_impl<E, tuple<Types...>> {
  using type = ind_expr<E>;
};

template <typename... Types>
struct is_indexed_impl<tuple<Types...>> : std::true_type {};

template <typename... Types> class tuple_wrapper {
public:
   tuple_wrapper(Types const &...values) : data{values...} {}

  template <typename E>
   tuple_wrapper(ind_expr<E> const &e) : data{e} {}

  template <size_t I>  decltype(auto) get() {
    return data.template get<I>();
  }

  template <size_t I>  decltype(auto) get() const {
    return data.template get<I>();
  }

  static constexpr size_t num_types = tuple<Types...>::num_types;

  template <size_t I>
  using ith_type = typename tuple<Types...>::template ith_type<I>;

private:
  tuple<Types...> data;
};
} // namespace nitro
