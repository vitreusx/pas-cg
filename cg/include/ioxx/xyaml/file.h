#pragma once
#include "node.h"

namespace ioxx::xyaml {
class file {
public:
  file() = default;

  std::string const &fetch();

public:
  std::optional<std::string> source;
  std::optional<std::filesystem::path> rel_path, abs_path;
};

template <> struct xyaml_conv<file> {
  void load(node const &from, file &to) const;
  void save(node &to, file const &from) const;
};

} // namespace ioxx::xyaml