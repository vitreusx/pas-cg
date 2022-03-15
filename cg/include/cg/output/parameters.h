#pragma once
#include <cg/utils/quantity.h>
#include <filesystem>
#include <ioxx/ioxx.h>

namespace cg::out {
struct parameters {
public:
  bool enabled;
  quantity stats_period, file_period;
  std::filesystem::path output_dir;

  void load(ioxx::xyaml::node const &from);
};
} // namespace cg::out