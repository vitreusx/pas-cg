#pragma once
#include "state.h"
#include "thread.h"
#include <vector>

namespace cg::simul {
class program {
public:
  static void main(int argc, char **argv);

private:
  static void perform_simulation(parameters const &params,
                                 ioxx::xyaml::node const &raw_params);
  static void check_determinism(parameters const &params,
                                ioxx::xyaml::node const &raw_params);
  static void run_from_checkpoint(std::filesystem::path const &ckpt_path);
};
} // namespace cg::simul