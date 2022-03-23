#pragma once
#include "../expr.h"
#include "../indexed.h"
#include "../macros.h"
#include "../tuple.h"

namespace nitro {
template <typename T> class set_node;

template <typename E> struct set_node_expr : ind_expr<E> {
  EXPR_BODY(item, is_vacant)

  decltype(auto) has_item() { return !is_vacant(); }
  decltype(auto) has_item() const { return !is_vacant(); }

  void remove() { is_vacant() = true; }
};

template <typename T> struct is_indexed_impl<set_node<T>> : std::true_type {};

template <typename E, typename T> struct expr_impl<E, set_node<T>> {
  using type = set_node_expr<E>;
};

template <typename E> struct set_node_auto_expr : set_node_expr<E> {
  AUTO_EXPR_BODY(item, is_vacant)
};

template <typename E, typename T> struct auto_expr_impl<E, set_node<T>> {
  using type = set_node_auto_expr<E>;
};

template <typename T> using set_node_base = nitro::tuple_wrapper<T, bool>;

template <typename T>
class set_node : public set_node_auto_expr<set_node<T>>,
                 public set_node_base<T> {
public:
  using Base = set_node_base<T>;
  using Base::Base;
  using Base::get;
};
} // namespace nitro