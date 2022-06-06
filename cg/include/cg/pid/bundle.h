#pragma once
#include <cg/types/amp.h>

namespace cg::pid {
template <typename E> struct bundle_expr { EXPR(i1, i2, orig_dist, type) };

class bundle : public bundle_expr<bundle> {
public:
  INST(bundle, FIELD(int, i1), FIELD(int, i2), FIELD(real, orig_dist),
       FIELD(int16_t, type))
};
} // namespace cg::pid
