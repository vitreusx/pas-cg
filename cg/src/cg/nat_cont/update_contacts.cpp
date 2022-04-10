#include <cg/base_forces/lj.h>
#include <cg/nat_cont/update_contacts.h>

namespace cg::nat_cont {

void update_contacts::operator()() const {
  contacts->clear();

  for (int idx = 0; idx < all_contacts.size(); ++idx) {
    auto nat_cont = all_contacts[idx];
    auto idx1 = nat_cont.i1(), idx2 = nat_cont.i2();
    auto nat_dist = nat_cont.nat_dist();

    auto cutoff = lj::compute_cutoff(nat_dist);

    auto r1 = r[idx1], r2 = r[idx2];
    auto cur_dist = norm(simul_box->r_uv(r1, r2));
    if (cur_dist < cutoff + nl->orig_pad) {
      contacts->push_back(nat_cont);
    }
  }
}
} // namespace cg::nat_cont