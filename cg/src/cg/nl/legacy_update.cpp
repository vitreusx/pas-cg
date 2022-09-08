#include <cg/nl/legacy_update.h>

namespace cg::nl {
template <typename E1, typename E2>
static bool operator<(nl::pair_expr<E1> const &left,
                      nl::pair_expr<E2> const &right) {
  return std::make_pair(left.i1(), left.i2()) <
         std::make_pair(right.i1(), right.i2());
}

void legacy_update::omp_async() const {
  auto req_r = cutoff + pad;

#pragma omp master
  nl_data->pairs.clear();

#pragma omp barrier

  int num_pairs = num_particles * num_particles;

#pragma omp for schedule(static) nowait
  for (int flat = 0; flat < num_pairs; ++flat) {
    int i1 = flat / num_particles, i2 = flat % num_particles;
    if (i1 >= i2)
      continue;

    auto r1 = r[i1];
    auto chain1 = chain_idx[i1], seq1 = seq_idx[i1];

    auto r2 = r[i2];
    auto chain2 = chain_idx[i2], seq2 = seq_idx[i2];

    auto diff = seq1 > seq2 ? seq1 - seq2 : seq2 - seq1;
    if (chain1 == chain2 && diff < 3)
      continue;

    auto orig_dist = norm(simul_box->wrap<vec3r>(r1, r2));
    if (orig_dist < req_r)
      nl_data->pairs.emplace_back(i1, i2, orig_dist, false);
  }
#pragma omp barrier

#pragma omp master
  {
    std::sort(nl_data->pairs.begin(), nl_data->pairs.end());

    nl_data->orig_pad = pad;
    for (int idx = 0; idx < r.size(); ++idx)
      nl_data->orig_r[idx] = r[idx];
    nl_data->orig_pbc = *simul_box;
    *invalid = false;
  }
}

} // namespace cg::nl