#pragma once
#include "macros.h"

namespace nitro {
template <typename E>
struct set_node_expr {
  EXPR(item, is_vacant);

  decltype(auto) has_item() const {
    return !is_vacant();
  }

  void remove() {
    is_vacant() = true;
  }
};

template <typename T>
struct set_node : set_node_expr<set_node<T>> {
  INST(set_node, FIELD(T, item), FIELD(bool, is_vacant));
};
} // namespace nitro