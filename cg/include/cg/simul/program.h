#pragma once
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
  void regular_main();
  void determinism_main();
  void setup_omp();
};
} // namespace cg::simul