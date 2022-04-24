#include <cg/files/sl4.h>

namespace ioxx::sl4 {
std::ostream &operator<<(std::ostream &os, element const &el) {
  el.write(os);
  return os;
}

void div::write(std::ostream &os) const {
  bool first = true;
  for (auto const &child : children) {
    if (!first)
      os << '\n';
    os << *child;
    first = false;
  }
}

element *div::find(const std::string &id) {
  if (auto iter = named_children.find(id); iter != named_children.end())
    return iter->second;
  else
    return nullptr;
}

element &div::operator[](int idx) {
  if (idx >= children.size())
    throw std::runtime_error("index out of bounds");
  else
    return *children[idx];
}

element &div::operator[](const std::string &id) {
  auto *el = find(id);
  if (el != nullptr)
    return *el;
  else
    throw std::runtime_error("children with id \"" + id + "\" not found");
}

table::table(ioxx::table::table tab) : tab{tab} {}

void table::write(std::ostream &os) const {
  ioxx::table::sl4_parser parser;
  parser.write(os, tab);
}

ioxx::table::table *table::operator->() { return &tab; }

comment::comment(std::string text) {
  std::for_each(text.begin(), text.end(),
                [](char &c) { c = (char)std::toupper(c); });
  this->text = "#" + std::move(text);
}

void comment::write(std::ostream &os) const { os << text; }

raw::raw(std::string text) : text{std::move(text)} {}

void raw::write(std::ostream &os) const { os << text; }
} // namespace ioxx::sl4