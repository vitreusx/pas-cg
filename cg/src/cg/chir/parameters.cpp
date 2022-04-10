#include <cg/chir/parameters.h>
#include <cg/utils/ioxx_interop.h>
namespace cg::chir {

void parameters::load(ioxx::xyaml::node const &p) {
  enabled = p["enabled"].as<bool>();
  e_chi = p["e_chi"].as<quantity>();
}
} // namespace cg::chir