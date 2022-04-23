#include <cg/files/table/parser.h>
#include <sstream>

namespace ioxx::table {
table parser::load(const std::string &source) const {
  std::stringstream ss{};
  ss << source;
  return read(ss);
}

std::string parser::dump(table const &tab) const {
  std::stringstream ss{};
  write(ss, tab);
  return ss.str();
}
} // namespace ioxx::table