#pragma once
#include <cg/simul/state.h>
#include <fstream>
#include <memory>

namespace cg::out {
class print_raw_data {
public:
  std::shared_ptr<std::ofstream> data_file;
  simul::state const *st;

public:
  void operator()() const;
};
} // namespace cg::out