#include <cg/tether/parameters.h>
#include <cg/utils/ioxx_interop.h>
namespace cg::tether {

void parameters::load(ioxx::xyaml::node const &p) {
  enabled = p["enabled"].as<bool>();
  H1 = p["H1"].as<quantity>();
  H2 = p["H2"].as<quantity>();
  def_length = p["default bond length"].as<quantity>();
}
} // namespace cg::tether