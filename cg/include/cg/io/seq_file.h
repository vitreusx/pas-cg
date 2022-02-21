#pragma once
#include "cg/state/model.h"
#include <filesystem>
#include <ioxx/xyaml.h>

namespace cg {
class seq_file {
public:
  seq_file() = default;
  explicit seq_file(std::filesystem::path const &seq_file_path);
  model xmd_model;

  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg