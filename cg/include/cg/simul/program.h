#pragma once
#include "state.h"
#include "thread.h"
#include <vector>

namespace cg::simul {
class program {
public:
  void main(int argc, char **argv);

private:
  state st;

  void parse_args(int argc, char **argv);
  void setup_omp();
};
} // namespace cg::simul