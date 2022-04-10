#include <cg/utils/ioxx_interop.h>
namespace ioxx::xyaml {
void xyaml_conv<cg::quantity>::load(const node &from, cg::quantity &to) const {
  std::string repr;
  from >> repr;
  to = cg::quantity(repr);
}

void xyaml_conv<cg::quantity>::save(node &to, const cg::quantity &from) const {
  auto repr = from.repr();
  to << repr;
}
} // namespace ioxx::xyaml
