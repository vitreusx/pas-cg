#include "xyaml/file.h"
#include <fstream>
#include <stdexcept>
using namespace ioxx::xyaml;

void xyaml_conv<file>::load(const node &from, file &to) const {
  if (from.IsScalar()) {
    to.source = from.as<std::string>();
    to.rel_path = to.abs_path = std::nullopt;
  } else if (from["(at path)"]) {
    to.rel_path = from["(at path)"].as<std::string>();
    to.abs_path = from["(at path)"].abs_path(to.rel_path.value());
  } else {
    throw std::runtime_error("invalid schema");
  }
}

void xyaml_conv<file>::save(node &to, const file &from) const {
  if (from.rel_path.has_value()) {
    to["(at path)"] = from.rel_path.value().string();
    auto path = to.abs_path(from.rel_path.value());

    if (from.source.has_value()) {
      std::filesystem::create_directories(path.parent_path());
      std::ofstream file(path);
      file << from.source.value();
    }
  } else {
    to = from.source.value();
  }
}

std::string const &file::fetch() {
  if (!source.has_value()) {
    std::ifstream source_file(abs_path.value());
    std::stringstream file_ss;
    file_ss << source_file.rdbuf();
    source = file_ss.str();
  }
  return source.value();
}