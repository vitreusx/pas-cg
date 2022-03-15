#include "output/make_report.h"
#include "utils/text.h"
namespace cg::out {

void make_report::operator()() const {
  using namespace ioxx::xyaml;

  bool gen_report = false;

  if (state->simul_first_time) {
    std::filesystem::remove_all(state->out_dir);
    auto simul_path = state->out_dir / "report.yml";
    state->for_simul = node::new_file(simul_path);
  }

  auto traj_dir = format("traj-%d", *state->traj_idx);
  if (state->traj_first_time) {
    auto traj_path = state->out_dir / traj_dir / "report.yml";
    state->for_traj = node::new_file(traj_path);
    gen_report = true;
  }

  gen_report |= (*t - *last_t >= period);

  if (gen_report) {
    auto step_dir = format("step-%d", state->step_idx);
    auto step_path = state->out_dir / traj_dir / step_dir / "report.yml";
    state->for_step = node::new_file(step_path);

    for (auto const *hook : *hooks)
      hook->report_to(*state);

    state->for_simul.save();
    state->for_traj.save();
    state->for_step.save();

    ++state->step_idx;
    *last_t = *t;
  }
  
  state->traj_first_time = false;
  state->simul_first_time = false;
}
} // namespace cg::out