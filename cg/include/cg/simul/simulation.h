#pragma once
#include "dynamics.h"
#include "kernels.h"
#include "parameters.h"
#include "thread_state.h"
#include <cg/amino/compiled.h>
#include <cg/input/model.h>
#include <cg/utils/random.h>
#include <fstream>
#include <iostream>
#include <mutex>
#include <omp.h>
#include <thread>

namespace cg::simul {
class simulation {
private:
  static void print_help(char **argv);

  std::vector<std::string> param_paths;
  void parse_args(int argc, char **argv);

  parameters params;
  void load_parameters();

  rand_gen gen;
  void general_setup();

  input::model orig_model, model;
  void load_model();

  input::model::res_map_t res_map;
  int num_res;
  nitro::vector<vec3r> orig_r, r;
  nitro::vector<amino_acid> atype;
  compiled_aa_data comp_aa_data;
  cg::box<real> box;
  nitro::vector<int> prev, next, chain_idx, seq_idx;
  void compile_model();

  real t, V;
  dynamics dyn;
  void setup_dyn();

  kernels ker;

  real report_last_t;
  out::report_state report;
  std::vector<out::hook const *> hooks;
  void setup_output();

  nitro::vector<real> mass_inv, mass_rsqrt;
  nitro::vector<vec3r> v;
  nitro::vector<vec3sr> y0, y1, y2, y3, y4, y5;
  solver_real true_t;
  void setup_langevin();

  bool pbar_first_time;
  pbar::render::time_point_t pbar_start_clock, pbar_last_clock;
  void setup_pbar();

  nl::data nl;
  bool nl_required, nl_invalid, verify_first_time;
  real max_cutoff;
  nitro::vector<int> res_cell_idx, reordered_idx, num_res_in_cell, cell_offset;
  nitro::vector<nl::pair> all_pairs;
  void setup_nl();

  nitro::vector<chir::chiral_quad> chir_quads;
  void setup_chir();

  nitro::vector<tether::pair> tether_pairs;
  void setup_tether();

  nitro::vector<nat_ang::nat_ang> native_angles;
  void setup_nat_ang();

  nitro::vector<heur_ang::heur_ang> heur_angles;
  void setup_heur_ang();

  nitro::vector<nat_dih> native_dihedrals;
  void setup_nat_dih();

  nitro::vector<heur_dih::heur_dih> heur_dihedrals;
  void setup_heur_dih();

  nitro::vector<pauli::pair> pauli_pairs;
  void setup_pauli();

  nitro::vector<nat_cont::nat_cont> all_native_contacts, cur_native_contacts;
  nitro::vector<nl::exclusion> nat_cont_excl;
  void setup_nat_cont();

  nitro::vector<dh::pair> dh_pairs;
  void setup_dh();

  nitro::set<qa::free_pair> qa_free_pairs;
  nitro::vector<qa::candidate> qa_candidates;
  nitro::set<qa::contact> qa_contacts;
  nitro::vector<sync_data> sync_values;
  nitro::vector<vec3r> n, h;
  nitro::vector<qa::cys_neigh> qa_cys_neigh;
  bool ss_spec_crit;
  nitro::vector<int> neigh_count, cys_indices, qa_removed;
  nitro::vector<bool> part_of_ssbond;
  int num_qa_contacts;
  void setup_qa();

  nitro::vector<pid::bundle> pid_bundles;
  nitro::vector<sink_lj> ss_ljs;
  void setup_pid();

  nitro::vector<vafm::tip> vafm_tips;
  void setup_vafm();

  nitro::vector<fafm::tip> fafm_tips;
  void setup_fafm();

  void setup(int argc, char **argv);

  class thread {
  public:
    thread() = default;
    explicit thread(simulation *simul);
    void main_async();

  private:
    dynamics dyn;
    rand_gen gen;
    kernels ker;
    parameters params;
    simulation *simul;

    void async_part();
    void sync_part(bool pre_loop = false);
    void on_nl_invalidation();
  };

public:
  int operator()(int argc, char **argv);
};
} // namespace cg::simul