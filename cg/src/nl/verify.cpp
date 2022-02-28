#include "nl/verify.h"
using namespace cg::nl;

void verify::operator()() const {
  if (*first_time) {
    *first_time = false;
    *invalid = true;
  } else {
    real max_r_disp = 0.0f;
    for (int idx = 0; idx < num_particles; ++idx) {
      auto r_ = r[idx], orig_r_ = nl_data->orig_r.at(idx);
      auto dr = simul_box->r_uv(r_, orig_r_);
      max_r_disp = max(max_r_disp, norm(dr));
    }

    auto box_disp = l1_norm(simul_box->cell - nl_data->orig_box.cell);

    real total_disp = 2.0f * max_r_disp + box_disp;
    *invalid = (total_disp > nl_data->orig_pad);
  }
}