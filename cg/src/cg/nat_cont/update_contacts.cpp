#include <cg/base_forces/lj.h>
#include <cg/nat_cont/update_contacts.h>
#include <iostream>

namespace cg::nat_cont {

void update_contacts::operator()() const {
  contacts->clear();

  auto nat_cont_iter = all_contacts.begin();

  for (decltype(auto) nl_pair : nl->pairs) {
    if (nl_pair.taken())
      continue;

    auto nl_pair_ = std::make_pair(nl_pair.i1(), nl_pair.i2());

    while (nat_cont_iter != all_contacts.end()) {
      auto nat_cont = *nat_cont_iter;
      auto nat_pair = std::make_pair(nat_cont.i1(), nat_cont.i2());

      if (nat_pair < nl_pair_) {
        ++nat_cont_iter;
      } else {
        if (nat_pair == nl_pair_) {
          auto r1 = r[nl_pair.i1()], r2 = r[nl_pair.i2()];
          auto cur_dist = norm(simul_box->wrap<vec3r>(r1, r2));
          if (cur_dist < cutoff + nl->orig_pad) {
            contacts->push_back(nat_cont);
            nl_pair.taken() = true;
          }
        }
        break;
      }
    }
  }
}

void update_contacts::mark_as_taken(nl::data *nl_) const {
  auto nat_cont_iter = all_contacts.begin();

  for (decltype(auto) nl_pair : nl_->pairs) {
    if (nl_pair.taken())
      continue;

    auto nl_pair_ = std::make_pair(nl_pair.i1(), nl_pair.i2());

    while (nat_cont_iter != all_contacts.end()) {
      auto nat_cont = *nat_cont_iter;
      auto nat_pair = std::make_pair(nat_cont.i1(), nat_cont.i2());

      if (nat_pair < nl_pair_) {
        ++nat_cont_iter;
      } else {
        if (nat_pair == nl_pair_) {
          auto r1 = r[nl_pair.i1()], r2 = r[nl_pair.i2()];
          auto cur_dist = norm(simul_box->wrap<vec3r>(r1, r2));
          if (cur_dist < cutoff + nl_->orig_pad) {
            nl_pair.taken() = true;
          }
        }
        break;
      }
    }
  }
}
} // namespace cg::nat_cont