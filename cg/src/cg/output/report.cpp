#include <cg/output/report.h>
#include <fstream>
#include <iostream>
namespace cg::out {
auto with_ext(std::filesystem::path const &prefix, std::string const &ext) {
  return prefix.parent_path() / (prefix.filename().string() + ext);
}

std::ofstream open_file(std::filesystem::path const &path) {
  if (path.has_parent_path()) {
    auto par = path.parent_path();
    if (!exists(par))
      std::filesystem::create_directories(par);
  }

  return std::ofstream(path);
}

void make_report::operator()() const {
  if (*t >= rep->stats_last + stats_every) {
    // emit stats
    rep->stats_last = *t;
  }

  if (*t >= rep->struct_last + struct_every) {
    auto xmd_model = *model;
    for (int idx = 0; idx < xmd_model.residues.size(); ++idx) {
      auto *ref_res = model->residues[idx].get();
      auto res_idx = res_map->at(ref_res);
      xmd_model.residues[idx]->pos = r->at(res_idx);
    }

    auto &pdb_model = rep->full_pdb.find_or_add_model(rep->model_serial);
    pdb_model = pdb_file(xmd_model).primary_model();
    pdb_model.model_serial = rep->model_serial;
    ++rep->model_serial;

    auto pdb_path = with_ext(prefix, ".pdb");
    auto pdb_of = open_file(pdb_path);
    pdb_of << rep->full_pdb;

    rep->struct_last = *t;
  }
}
} // namespace cg::out