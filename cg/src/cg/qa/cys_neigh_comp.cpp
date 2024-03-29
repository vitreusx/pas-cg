#include <cg/qa/cys_neigh_comp.h>
namespace cg::qa {

void update_cys_neigh::operator()() const {
  neigh->clear();

  for (decltype(auto) nl_pair : nl->pairs) {
    auto i1 = nl_pair.i1(), i2 = nl_pair.i2();
    auto r1 = r[i1], r2 = r[i2];

    if (norm(simul_box->wrap<vec3r>(r1, r2)) < neigh_radius + nl->orig_pad) {
      neigh->emplace_back(i1, i2);
    }
  }
}

void reset_cys_neigh::operator()() const {
  for (int idx = 0; idx < cys_indices.size(); ++idx)
    neigh_count[cys_indices[idx]] = 0;
}

void count_cys_neigh::iter(int idx) const {
  auto pair = neigh->at(idx);
  auto r_cys = r[pair.cys_idx()], r_nei = r[pair.neigh_idx()];
  if (norm(simul_box->wrap<vec3r>(r_cys, r_nei)) < neigh_radius)
    ++neigh_count[pair.cys_idx()];
}

void count_cys_neigh::operator()() const {
  for_slice(0, total_size());
}

void count_cys_neigh::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < neigh->size(); ++idx)
    iter(idx);
}

void count_cys_neigh::for_slice(int from, int to) const {
  for (int idx = from; idx < to; ++idx)
    iter(idx);
}

int count_cys_neigh::total_size() const {
  return neigh->size();
}

} // namespace cg::qa