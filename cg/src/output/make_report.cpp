#include "output/make_report.h"
namespace cg::out {

void make_report::operator()() const {

  bool gen_report = false;

  if (state->simul_first_time)
    state->simul_init();

  if (state->traj_first_time) {
    state->traj_init();
    gen_report = true;
  }

  gen_report |= (*t - *last_t >= period);

  if (gen_report) {
    state->step_init();
    for (auto const *hook : *hooks)
      hook->report_to(*state);
    state->step_finish();

    *last_t = *t;
  }

  state->traj_first_time = false;
  state->simul_first_time = false;
}
} // namespace cg::out