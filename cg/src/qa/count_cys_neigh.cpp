#include "qa/count_cys_neigh.h"
using namespace cg::qa;

void update_cys_neigh::operator()() const {
  neigh->clear();

  for (int pair_idx = 0; pair_idx < nl->non_native.size(); ++pair_idx) {
    auto nl_pair = nl->non_native[pair_idx];
    auto i1 = nl_pair.i1(), i2 = nl_pair.i2();
    auto r1 = r[i1], r2 = r[i2];

    if (norm(simul_box->r_uv(r1, r2)) < neigh_radius + nl->orig_pad) {
      neigh->emplace_back(i1, i2);
    }
  }
}

void count_cys_neigh::operator()() const {
  for (int idx = 0; idx < cys_indices.size(); ++idx) {
    neigh_count[cys_indices[idx]] = 0;
  }

  for (int idx = 0; idx < neigh->size(); ++idx) {
    auto pair = neigh->at(idx);
    auto r_cys = r[pair.cys_idx()], r_nei = r[pair.neigh_idx()];
    if (norm(simul_box->r_uv(r_cys, r_nei)) < neigh_radius)
      ++neigh_count[pair.cys_idx()];
  }
}