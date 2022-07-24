#pragma once
#include "parameters.h"
#include <cg/input/pdb_file.h>
#include <cg/types/amp.h>
#include <fstream>
#include <limits>
#include <memory>

namespace cg::out {
class legacy_report {
public:
  std::shared_ptr<std::fstream> out_file, map_file, pdb_file;
  real stats_last, struct_last;
};
} // namespace cg::out