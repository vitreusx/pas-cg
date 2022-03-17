#include "output/make_report.h"
namespace cg::out {

void make_report::operator()() const {
  if (state->simul_first_time)
    state->setup_simul_report();

  if (state->traj_first_time)
    state->setup_traj_report();

  state->report_stats =
      state->traj_first_time || (*t - state->last_stats_t >= stats_period);
  state->report_files =
      state->traj_first_time || (*t - state->last_files_t >= file_period);

  if (state->report_stats || state->report_files) {
    state->setup_step_report();
    for (auto const *hook : *hooks)
      hook->report_to(*state);
    state->finalize_step_report();

    if (state->report_stats)
      state->last_stats_t = *t;
    if (state->report_files)
      state->last_files_t = *t;
  }
}
} // namespace cg::out