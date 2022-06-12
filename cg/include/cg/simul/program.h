#pragma once
#include "legacy.h"
#include "state.h"
#include "thread.h"
#include <vector>

namespace cg::simul {
class program {
public:
  void main(int argc, char **argv);

private:
  parameters params;

  void parse_args(int argc, char **argv);
  void perform_simulation();
  void check_determinism();
  void run_legacy_version();
  void setup_omp();
};
} // namespace cg::simul