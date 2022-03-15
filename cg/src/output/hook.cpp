#include "output/hook.h"

namespace cg::out {
void report_state::on_new_trajectory() {
  traj_first_time = true;
  step_idx = 0;
}
} // namespace cg::out