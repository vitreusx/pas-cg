#pragma once
#include "input/model.h"
#include <filesystem>
#include <ioxx/xyaml.h>

namespace cg {
class seq_file {
public:
  seq_file() = default;
  input::model model;

  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg