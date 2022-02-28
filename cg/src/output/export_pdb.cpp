#include "output/export_pdb.h"
#include "input/pdb_file.h"
#include <fstream>
using namespace cg::out;

void export_pdb::report_to(report_state &report) const {
  using namespace ioxx::xyaml;

  auto cur_model = *ref_model;

  for (int idx = 0; idx < (int)cur_model.residues.size(); ++idx) {
    auto *ref_res = ref_model->residues[idx].get();
    auto res_idx = res_map->at(ref_res);
    cur_model.residues[idx]->pos = r[res_idx];
  }

  auto pdb_repr = pdb_file(cur_model);
  std::stringstream pdb_ss;
  pdb_ss << pdb_repr.emit_model(report.ord + 1);

  if (report.first_time) {
    file gen_model_file;
    gen_model_file.source = "";
    gen_model_file.rel_path = "model.pdb";
    report.general["model"] = gen_model_file;
  }

  auto gen_model_path = report.general.abs_path("model.pdb");
  std::ofstream gen_model_file(gen_model_path, std::ios::app);
  if (!report.first_time)
    gen_model_file << '\n';
  gen_model_file << pdb_ss.str();

  file coords_file;
  coords_file.source = pdb_ss.str();
  coords_file.rel_path = "model.pdb";
  report.current["model"] = coords_file;
}