#include <cg/files/xyaml/file.h>
#include <cg/files/xyaml/node.h>

namespace ioxx::xyaml {
YAML::Node node::flatten() const {
  if (IsMap() && this->operator[]("(at path)")) {
    return YAML::Node(as<file>().fetch());
  } else if (IsMap()) {
    YAML::Node flat;
    for (auto iter = begin(); iter != end(); ++iter) {
      flat[iter->first] = this->operator[](iter->first).flatten();
    }
    if (children) {
      for (auto const &[key, val] : *children) {
        flat[key] = val.flatten();
      }
    }
    return flat;
  } else {
    return YAML::Clone(*this);
  }
}
} // namespace ioxx::xyaml