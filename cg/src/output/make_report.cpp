#include "output/make_report.h"
namespace cg::out {

void make_report::operator()() const {
  if (state->simul_first_time)
    state->simul_init();

  if (state->traj_first_time)
    state->traj_init();

  state->report_stats =
      state->traj_first_time || (*t - state->last_stats_t >= stats_period);
  state->report_files =
      state->traj_first_time || (*t - state->last_files_t >= file_period);

  if (state->report_stats || state->report_files) {
    state->step_init();

    for (auto const *hook : *hooks)
      hook->report_to(*state);
    state->step_finish();

    if (state->report_stats)
      state->last_stats_t = *t;
    if (state->report_files)
      state->last_files_t = *t;
  }

  state->traj_first_time = false;
  state->simul_first_time = false;
}
} // namespace cg::out