#include "features/nat_cont/update_contacts.h"
#include "base_forces/lj.h"
using namespace cg::nat_cont;

void update_contacts::operator()() const {
  contacts->clear();

  for (int idx = 0; idx < all_contacts->size(); ++idx) {
    auto nat_cont = all_contacts->at(idx);
    auto idx1 = nat_cont.i1(), idx2 = nat_cont.i2();
    auto nat_dist = nat_cont.nat_dist();

    auto cutoff = lj::cutoff(nat_dist);

    auto r1 = r->at(idx1), r2 = r->at(idx2);
    if (norm(box->r_uv(r1, r2)) < cutoff + nl->orig_pad) {
      contacts->push_back(nat_cont);
    }
  }
}