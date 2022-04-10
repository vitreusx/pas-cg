#include <cg/nat_dih/complex/parameters.h>
#include <cg/utils/ioxx_interop.h>
namespace cg::cnd {

void parameters::load(ioxx::xyaml::node const &p) {
  enabled = p["enabled"].as<bool>();
  CDA = p["CDA"].as<quantity>();
  CDB = p["CDB"].as<quantity>();
}
} // namespace cg::cnd