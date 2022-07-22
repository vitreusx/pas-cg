#pragma once
#include "parameters.h"
#include <cg/input/pdb_file.h>
#include <cg/types/amp.h>
#include <limits>

namespace cg::out {
struct snapshot {
  real t;
  vect::vector<vec3r> r;
  sbox::pbc<real> model_box;
};

class report {
public:
  real stats_last, struct_last;

  std::vector<snapshot> snapshots;
  ioxx::sl4::div out_root, map_root;

  void traj_init(int traj_idx);
};
} // namespace cg::out