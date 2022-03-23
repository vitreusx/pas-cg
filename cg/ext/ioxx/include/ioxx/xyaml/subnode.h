#pragma once
#include "node.h"

namespace ioxx::xyaml {
class subnode : public node {
public:
  using node::node;
  using node::operator=;

  subnode(node const &base);
  subnode &operator=(node const &base);
};

template <> struct xyaml_conv<subnode> {
  void load(node const &from, subnode &to) const;
  void save(node &to, subnode const &from) const;
};
} // namespace ioxx::xyaml