#pragma once
#include <cg/input/model.h>
#include <filesystem>
#include <ioxx/ioxx.h>

namespace cg {
class seq_file {
public:
  seq_file() = default;
  input::model model;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg