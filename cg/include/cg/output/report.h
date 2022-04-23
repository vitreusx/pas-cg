#pragma once
#include "parameters.h"
#include <cg/input/pdb_file.h>
#include <cg/types/amp.h>
#include <limits>

namespace cg::out {
class report {
public:
  real stats_last, struct_last;
  int model_serial;
  pdb_file full_pdb;
};

class make_report {
public:
  std::filesystem::path prefix;
  real stats_every, struct_every;

public:
  report *rep;
  real const *t;
  input::model const *model;
  input::model::res_map_t const *res_map;
  nitro::vector<vec3r> const *r;

public:
  void operator()() const;
};
} // namespace cg::out