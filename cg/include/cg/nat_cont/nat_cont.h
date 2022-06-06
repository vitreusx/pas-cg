#pragma once
#include "type.h"
#include <cg/types/amp.h>
#include <cg/vect/vect.h>

namespace cg::nat_cont {
template <typename E> struct nat_cont_expr {
  EXPR(i1, i2, nat_dist, type, formed, formation_t, all_cont_idx)

  decltype(auto) is_ssbond() const {
    return type() == type::SSBOND;
  }
};

class nat_cont : public nat_cont_expr<nat_cont> {
public:
  INST(nat_cont, FIELD(int, i1), FIELD(int, i2), FIELD(real, nat_dist),
       FIELD(cg::nat_cont::type, type), FIELD(bool, formed),
       FIELD(real, formation_t), FIELD(int, all_cont_idx))
};
} // namespace cg::nat_cont
