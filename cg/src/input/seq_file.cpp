#include "input/seq_file.h"
#include "map_file.h"
#include <ioxx/xyaml.h>
#include <regex>
#include <set>

namespace cg {
void seq_file::connect(ioxx::xyaml_proxy &proxy) {
  using namespace ioxx;

  auto file_node = xyaml_embedded();
  proxy &file_node;
  auto file_proxy = proxy(file_node.node);

  if (proxy.loading()) {
    int chain_idx = 0;

    for (auto const &chain_node : file_proxy["model"]["chains"]) {
      auto res_codes = chain_node["seq"].as<std::string>();
      res_codes = std::regex_replace(res_codes, std::regex("\\s+"), "");

      auto &xmd_chain =
          model.chains.emplace_back(std::make_unique<input::model::chain>());
      xmd_chain->residues.resize(res_codes.size());
      xmd_chain->chain_idx = chain_idx++;

      for (size_t idx = 0; idx < res_codes.size(); ++idx) {
        auto res_code = res_codes[idx];

        auto &xmd_res = model.residues.emplace_back(
            std::make_unique<input::model::residue>());
        xmd_res->seq_idx = (int)idx;
        xmd_res->parent = &*xmd_chain;
        xmd_res->type = amino_acid(res_code);
        xmd_res->pos = {0.0, 0.0, 0.0};

        xmd_chain->residues[idx] = &*xmd_res;
      }

      std::set<std::tuple<input::model::residue *, input::model::residue *,
                          input::model::residue *>>
          taken_angles;
      std::set<std::tuple<input::model::residue *, input::model::residue *,
                          input::model::residue *, input::model::residue *>>
          taken_dihedrals;

      auto maps_node = chain_node["maps"];
      for (auto const &entry : maps_node) {
        auto xentry = xyaml_node::from_data(entry, file_proxy.location);
        auto xproxy = proxy(xentry);

        map_file mf;
        if (xproxy["__source"]) {
          xyaml_embedded mf_node;
          xproxy["__source"] & mf_node;
          mf = mf_node.node.as<map_file>();
          auto shift_val = xproxy["shift"].as<int>();
          mf.shift(shift_val);
        } else {
          xyaml_embedded mf_node;
          xproxy &mf_node;
          mf = mf_node.node.as<map_file>();
        }

        for (auto const &mf_cont : mf.contacts) {
          auto &xmd_cont = model.contacts.emplace_back();
          xmd_cont.res1 = xmd_chain->residues[mf_cont.i1];
          xmd_cont.res2 = xmd_chain->residues[mf_cont.i2];
          xmd_cont.length = mf_cont.length;
          xmd_cont.type = input::model::UNKNOWN;
        }

        for (auto const &mf_angle : mf.angles) {
          auto &xmd_angle = model.angles.emplace_back();
          xmd_angle.res1 = xmd_chain->residues[mf_angle.i1];
          xmd_angle.res2 = xmd_chain->residues[mf_angle.i2];
          xmd_angle.res3 = xmd_chain->residues[mf_angle.i3];
          xmd_angle.theta = mf_angle.theta;

          taken_angles.insert(
              std::make_tuple(xmd_angle.res1, xmd_angle.res2, xmd_angle.res3));
        }

        for (auto const &mf_dihedral : mf.dihedrals) {
          auto &xmd_dihedral = model.dihedrals.emplace_back();
          xmd_dihedral.res1 = xmd_chain->residues[mf_dihedral.i1];
          xmd_dihedral.res2 = xmd_chain->residues[mf_dihedral.i2];
          xmd_dihedral.res3 = xmd_chain->residues[mf_dihedral.i3];
          xmd_dihedral.res4 = xmd_chain->residues[mf_dihedral.i4];
          xmd_dihedral.phi = mf_dihedral.phi;

          taken_dihedrals.insert(
              std::make_tuple(xmd_dihedral.res1, xmd_dihedral.res2,
                              xmd_dihedral.res3, xmd_dihedral.res4));
        }
      }

      for (size_t idx = 0; idx + 1 < xmd_chain->residues.size(); ++idx) {
        auto *res1 = xmd_chain->residues[idx];
        auto *res2 = xmd_chain->residues[idx + 1];

        auto &xmd_tether = model.tethers.emplace_back();
        xmd_tether.res1 = res1;
        xmd_tether.res2 = res2;
        xmd_tether.length = std::nullopt;

        if (idx + 2 < xmd_chain->residues.size()) {
          auto *res3 = xmd_chain->residues[idx + 2];

          bool taken_angle =
              taken_angles.count(std::make_tuple(res1, res2, res3)) > 0;

          if (!taken_angle) {
            auto &xmd_angle = model.angles.emplace_back();
            xmd_angle.res1 = res1;
            xmd_angle.res2 = res2;
            xmd_angle.res3 = res3;
            xmd_angle.theta = std::nullopt;
          }

          if (idx + 3 < xmd_chain->residues.size()) {
            auto *res4 = xmd_chain->residues[idx + 3];

            bool taken_dihedral =
                taken_dihedrals.count(std::make_tuple(res1, res2, res3, res4)) >
                0;

            if (!taken_dihedral) {
              auto &xmd_dihedral = model.dihedrals.emplace_back();
              xmd_dihedral.res1 = res1;
              xmd_dihedral.res2 = res2;
              xmd_dihedral.res3 = res3;
              xmd_dihedral.res4 = res4;
              xmd_dihedral.phi = std::nullopt;
            }
          }
        }
      }
    }
  }
}
} // namespace cg