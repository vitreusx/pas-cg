#include "qa/finish_processing.h"
using namespace cg::qa;

void finish_processing::operator()() const {
  for (int idx = 0; idx < removed->size(); ++idx) {
    auto removed_idx = removed->at(idx);

    auto node = contacts->at(removed_idx);
    auto contact = node.item();
    auto i1 = contact.i1(), i2 = contact.i2();

    part_of_ssbond[i1] = false;
    part_of_ssbond[i2] = false;
    sync[i1] += contact.sync_diff1();
    sync[i2] += contact.sync_diff2();
    free_pairs->emplace(i1, i2);

    node.remove();
    --*num_contacts;
  }

  removed->clear();

  std::sort(candidates->begin(), candidates->end(),
            [](auto l, auto r) -> bool { return l.dist() < r.dist(); });

  for (int idx = 0; idx < candidates->size(); ++idx) {
    auto candidate = candidates->at(idx);

    auto i1 = candidate.i1(), i2 = candidate.i2();
    auto sync_diff1 = candidate.sync_diff1();
    auto sync_diff2 = candidate.sync_diff2();
    auto type = candidate.type();

    sync_data new_sync1 = sync[i1] - sync_diff1;
    sync_data new_sync2 = sync[i2] - sync_diff2;

    if (!new_sync1.is_valid() || !new_sync2.is_valid())
      return;

    static auto ss_type = contact_type::SIDE_SIDE(aa_code::CYS, aa_code::CYS);
    bool is_ssbond = (int16_t)type == (int16_t)ss_type;

    if (disulfide_special_criteria && is_ssbond) {
      if (part_of_ssbond[i1] || part_of_ssbond[i2])
        return;
    }

    contacts->emplace(i1, i2, type, FORMING_OR_FORMED, *t, sync_diff1,
                      sync_diff2);
    ++*num_contacts;

    sync[i1] = new_sync1;
    sync[i2] = new_sync2;

    free_pairs->remove(candidate.free_pair_idx());

    if (disulfide_special_criteria && is_ssbond) {
      part_of_ssbond[i1] = true;
      part_of_ssbond[i2] = true;
    }
  }

  candidates->clear();

  if (2 * *num_contacts < contacts->size())
    contacts->compress();
}