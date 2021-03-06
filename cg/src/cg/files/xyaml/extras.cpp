#include <cg/files/xyaml/extras.h>

namespace ioxx::xyaml {

void user_repr<std::filesystem::path>::load(const node &from,
                                            std::filesystem::path &to) const {
  to = from.as<std::string>();
}

void user_repr<std::filesystem::path>::save(
    node &to, const std::filesystem::path &from) const {
  to = from.string();
}

void user_repr<cg::quantity>::load(const node &from, cg::quantity &to) const {
  std::string repr;
  from >> repr;
  to = cg::quantity(repr);
}

void user_repr<cg::quantity>::save(node &to, const cg::quantity &from) const {
  auto repr = from.repr();
  to << repr;
}

} // namespace ioxx::xyaml