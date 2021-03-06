#include <cg/nl/legacy_update.h>
#include <cg/nl/order.h>
namespace cg::nl {

void legacy_update::operator()() const {
  auto req_r = cutoff + pad;

  nl_data->native.clear();
  nl_data->non_native.clear();

  int nat_cont_idx = 0;

  for (int i1 = 0; i1 < num_particles; ++i1) {
    auto r1 = r[i1];
    auto chain1 = chain_idx[i1], seq1 = seq_idx[i1];

    for (int i2 = i1 + 1; i2 < num_particles; ++i2) {
      auto r2 = r[i2];
      auto chain2 = chain_idx[i2], seq2 = seq_idx[i2];

      auto diff = seq1 > seq2 ? seq1 - seq2 : seq2 - seq1;
      if (chain1 == chain2 && diff < 3)
        continue;

      auto orig_dist = norm(simul_box->wrap(r1, r2));
      if (orig_dist < req_r) {
        bool non_native = true;
        auto cur_pair = pair(i1, i2, orig_dist);

        while (nat_cont_idx < all_nat_cont.size()) {
          auto cur_nat_cont = all_nat_cont[nat_cont_idx];

          if (cur_nat_cont < cur_pair) {
            ++nat_cont_idx;
          } else {
            if (!(cur_nat_cont > cur_pair)) {
              nl_data->native.push_back(cur_pair);
              non_native = false;
            }
            break;
          }
        }

        if (non_native) {
          nl_data->non_native.push_back(cur_pair);
        }
      }
    }
  }

  nl_data->orig_pad = pad;
  for (int idx = 0; idx < r.size(); ++idx)
    nl_data->orig_r[idx] = r[idx];
  nl_data->orig_pbc = *simul_box;
  *invalid = false;
}

void legacy_update::omp_async() const {
  auto req_r = cutoff + pad;

#pragma omp master
  {
    nl_data->native.clear();
    nl_data->non_native.clear();
  }
#pragma omp barrier

  int nat_cont_idx = 0;
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

    auto orig_dist = norm(simul_box->wrap(r1, r2));
    if (orig_dist < req_r) {
      bool non_native = true;
      auto cur_pair = pair(i1, i2, orig_dist);

      while (nat_cont_idx < all_nat_cont.size()) {
        auto cur_nat_cont = all_nat_cont[nat_cont_idx];

        if (cur_nat_cont < cur_pair) {
          ++nat_cont_idx;
        } else {
          if (!(cur_nat_cont > cur_pair)) {
#pragma omp critical
            nl_data->native.push_back(cur_pair);
            non_native = false;
          }
          break;
        }
      }

      if (non_native) {
#pragma omp critical
        nl_data->non_native.push_back(cur_pair);
      }
    }
  }
#pragma omp barrier

#pragma omp master
  {
    std::sort(nl_data->non_native.begin(), nl_data->non_native.end());
    std::sort(nl_data->native.begin(), nl_data->native.end());

    nl_data->orig_pad = pad;
    for (int idx = 0; idx < r.size(); ++idx)
      nl_data->orig_r[idx] = r[idx];
    nl_data->orig_pbc = *simul_box;
    *invalid = false;
  }
}

} // namespace cg::nl