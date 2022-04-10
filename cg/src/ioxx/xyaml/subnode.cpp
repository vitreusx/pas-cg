#include <fstream>
#include <ioxx/xyaml/subnode.h>

namespace ioxx::xyaml {

void xyaml_conv<subnode>::load(const node &from, subnode &to) const {
  if (from["(at path)"]) {
    auto rel_path = from["(at path)"].as<std::string>();
    auto abs_path = from.abs_path(rel_path);
    to = node::import(abs_path);
  } else {
    to = from;
  }
}

void xyaml_conv<subnode>::save(node &to, const subnode &from) const {
  if (to.loc != from.loc && from.loc.has_value()) {
    auto rel_path = to.rel_path(from.loc.value());
    to["(at path)"] = rel_path.string();

    std::ofstream subnode_file(from.loc.value());
    subnode_file << from;
  } else {
    to = from;
  }
}

subnode::subnode(const node &base) : node(base) {}

subnode &subnode::operator=(const node &base) {
  this->node::operator=(base);
  return *this;
}

} // namespace ioxx::xyaml