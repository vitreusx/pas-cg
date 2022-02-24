#pragma once
#include "dynamics.h"
#include "kernels.h"
#include "parameters.h"
#include <cg/amino/compiled.h>
#include <cg/input/model.h>
#include <cg/utils/random.h>
#include <fstream>
#include <iostream>
#include <omp.h>

namespace cg::simul {
class simulation {
private:
  static void print_help(char **argv);

  std::optional<std::string> param_path;
  void parse_args(int argc, char **argv);

  void setup();

  parameters params;
  void load_parameters();

  rand_gen gen;
  void general_setup();

  using res_map_t = std::unordered_map<input::model::model::residue *, int>;
  input::model model;
  void load_model();

  res_map_t res_map;
  int num_res;
  nitro::vector<vec3r> r;
  nitro::vector<amino_acid> atype;
  compiled_aa_data comp_aa_data;
  cg::box<real> box;
  nitro::vector<int> prev, next, chain_idx, seq_idx;
  void compile_model();

  real t, V;
  dynamics dyn;
  void setup_dyn();

  kernels ker;

  nitro::vector<real> mass_inv, mass_rsqrt;
  nitro::vector<vec3r> v;
  nitro::vector<vec3sr> y0, y1, y2, y3, y4, y5;
  solver_real true_t;
  void setup_langevin();

  bool pbar_first_time;
  pbar::render::time_point_t start_wall_time;
  real pbar_last_t;
  void setup_pbar();

  nl::data nl;
  bool nl_invalid, verify_first_time;
  real max_cutoff;
  void setup_nl();

  nitro::vector<chir::chiral_quad> chir_quads;
  void setup_chir();

  nitro::vector<tether::pair> tether_pairs;
  void setup_tether();

  nitro::vector<nat_ang::nat_ang> native_angles;
  void setup_nat_ang();

  nitro::vector<nat_dih> native_dihedrals;
  void setup_nat_dih();

  nitro::vector<pauli::pair> pauli_pairs;
  void setup_pauli();

  nitro::vector<nat_cont::nat_cont> all_native_contacts, cur_native_contacts;
  nitro::vector<nl::pair> native_contact_exclusions;
  void setup_nat_cont();

  nitro::vector<dh::pair> dh_pairs;
  void setup_dh();

  nitro::set<qa::free_pair> qa_free_pairs;
  nitro::vector<qa::candidate> qa_candidates;
  nitro::set<qa::contact> qa_contacts;
  nitro::vector<sync_data> sync_values;
  nitro::vector<vec3r> n, h;
  void setup_qa();

  nitro::vector<pid::bundle> pid_bundles;
  nitro::vector<sink_lj> ss_ljs;
  void setup_pid();

  nitro::vector<vafm::tip> vafm_tips;
  void setup_vafm();

  nitro::vector<fafm::tip> fafm_tips;
  void setup_fafm();

  thread_state state;

  void main_loop();
  void pre_loop_init();
  void on_nl_invalidation();
  void async_part();

public:
  int operator()(int argc, char **argv);
};
} // namespace cg::simul