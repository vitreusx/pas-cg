#pragma once
#include <cg/files/files.h>
#include <cg/utils/quantity.h>
#include <filesystem>

namespace cg::out {
struct parameters {
public:
  bool enabled;
  std::optional<quantity> stats_every, struct_every;
  std::filesystem::path prefix;

  void load(ioxx::xyaml::node const &from);
};
} // namespace cg::out