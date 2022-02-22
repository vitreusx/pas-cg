#include "features/qa/process_candidates.h"
using namespace cg::qa;

void process_candidates::operator()() const {
  for (int idx = 0; idx < candidates->size(); ++idx) {
    iter(candidates->at(idx));
  }
}

template <typename E>
void process_candidates::iter(candidate_expr<E> const &candidate) const {
  auto i1 = candidate.i1(), i2 = candidate.i2();
  auto sync_diff1 = candidate.sync_diff1();
  auto sync_diff2 = candidate.sync_diff2();
  auto type = candidate.type();

  sync_data new_sync1 = sync->at(i1) - sync_diff1;
  sync_data new_sync2 = sync->at(i2) - sync_diff2;

  if (new_sync1.is_valid() && new_sync2.is_valid()) {
    contacts->emplace_back(i1, i2, type, FORMING_OR_FORMED, *t, sync_diff1,
                           sync_diff2);
 
    sync->at(i1) = new_sync1;
    sync->at(i2) = new_sync2;

    free_pairs->remove(candidate.free_pair_idx());
  }
}