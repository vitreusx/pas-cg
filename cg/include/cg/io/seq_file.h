#pragma once
#include "cg/state/model.h"
#include <filesystem>
#include <ioxx/xyaml.h>

namespace cg {
class seq_file {
public:
  explicit seq_file(std::filesystem::path const &seq_file_path);
  bool connect(ioxx::xyaml_proxy &proxy);
  model const &to_model() const;

private:
  model xmd_model;
};
} // namespace cg