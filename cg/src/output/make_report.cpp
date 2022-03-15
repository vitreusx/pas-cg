#include "output/make_report.h"
namespace cg::out {

void make_report::operator()() const {
  using namespace ioxx::xyaml;
  if (state->first_time) {
    std::filesystem::remove_all(state->output_dir);
    auto gen_path = state->output_dir / "report.yml";
    state->general = node::new_file(gen_path);
  }

  if (state->first_time || *t - *last_t >= period) {
    auto cur_path =
        state->output_dir / std::to_string(state->ord) / "report.yml";
    state->current = node::new_file(cur_path);

    for (auto const *hook : *hooks)
      hook->report_to(*state);

    state->general.save();
    state->current.save();

    ++state->ord;
    *last_t = *t;
  }

  state->first_time = false;
}
} // namespace cg::out