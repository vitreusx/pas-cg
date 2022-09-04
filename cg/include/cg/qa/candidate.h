#pragma once
#include "contact_type.h"
#include <cg/amino/sync_data.h>
#include <cg/types/amp.h>

namespace cg::qa {
template <typename E>
struct candidate_expr {
  EXPR(i1, i2, orig_dist, dist, free_pair_idx, type, sync_diff1, sync_diff2)
};

class candidate : public candidate_expr<candidate> {
public:
  INST(candidate, FIELD(int, i1), FIELD(int, i2), FIELD(real, orig_dist),
       FIELD(real, dist), FIELD(int, free_pair_idx), FIELD(contact_type, type),
       FIELD(sync_data, sync_diff1), FIELD(sync_data, sync_diff2))
};
} // namespace cg::qa