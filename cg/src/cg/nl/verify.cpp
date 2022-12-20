#include <cg/nl/verify.h>

namespace cg::nl {
void verify::omp_async() const {
  if (*first_time) {
#pragma omp master
    {
      *first_time = false;
      *invalid = true;
    }
#pragma omp barrier
  } else {
    real max_r_disp = 0.0f;
#pragma omp for schedule(static) nowait
    for (auto idx : idxes) {
      auto r_ = r[idx], orig_r_ = nl_data->orig_r.at(idx);
      auto dr = simul_box->wrap<vec3r>(r_, orig_r_);
      max_r_disp = max(max_r_disp, norm(dr));
    }

    auto box_disp = l1_norm(simul_box->cell - nl_data->orig_pbc.cell);
#pragma omp single nowait
    *total_disp = (real)0.0f;

#pragma omp barrier

    real thread_disp = 2.0f * max_r_disp + box_disp;
#pragma omp critical
    {
      *total_disp = max(*total_disp, thread_disp);
      *invalid = (*total_disp > nl_data->orig_pad);
    }
  }
}

} // namespace cg::nl