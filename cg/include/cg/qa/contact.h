#pragma once
#include "contact_type.h"
#include <cg/amino/sync_data.h>
#include <cg/types/amp.h>

namespace cg::qa {
enum contact_status {
  FORMING_OR_FORMED,
  BREAKING
};

template <typename E> struct contact_expr {
  EXPR(i1, i2, orig_dist, type, status, ref_time, sync_diff1, sync_diff2)
};

class contact : public contact_expr<contact> {
public:
  INST(contact, FIELD(int, i1), FIELD(int, i2), FIELD(real, orig_dist),
       FIELD(contact_type, type), FIELD(contact_status, status),
       FIELD(real, ref_time), FIELD(sync_data, sync_diff1),
       FIELD(sync_data, sync_diff2));
};
} // namespace cg::qa