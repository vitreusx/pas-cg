#include "output/hook.h"
#include "utils/text.h"

namespace cg::out {
void report_data::simul_init() {
  using namespace ioxx::xyaml;
  std::filesystem::remove_all(out_dir);
  auto simul_path = out_dir / "report.yml";
  for_simul = node::new_file(simul_path);
}

void report_data::traj_init() {
  using namespace ioxx::xyaml;
  snap_idx = 0;
  traj_dir = out_dir / format("traj-%d", *traj_idx);
  auto traj_path = traj_dir / "report.yml";
  for_traj = node::new_file(traj_path);
}

void report_data::step_init() {
  using namespace ioxx::xyaml;
  snap_dir = traj_dir / format("snap-%d", snap_idx);
  auto step_path = snap_dir / "report.yml";
  for_snap = node::new_file(step_path);
  scalars.emplace_back();
}

void report_data::step_finish() {
  using namespace ioxx::xyaml;
  csv<> scalars_csv;
  scalars_csv.path = "scalars.csv";
  scalars_csv.data.header = ioxx::csv_header(columns);

  bool scalars_present = false;
  for (auto const &values : scalars) {
    auto &row = scalars_csv.data.emplace_back();
    for (auto const &[key, value] : values) {
      row[key] = value;
      scalars_present = true;
    }
  }

  if (scalars_present)
    for_traj["scalars"] = scalars_csv;

  for_simul.save();
  for_traj.save();
  for_snap.save();
  ++snap_idx;
}

} // namespace cg::out