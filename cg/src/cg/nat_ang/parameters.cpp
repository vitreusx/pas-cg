#include <cg/nat_ang/parameters.h>
#include <cg/utils/ioxx_interop.h>
namespace cg::nat_ang {

void parameters::load(ioxx::xyaml::node const &p) {
  enabled = p["enabled"].as<bool>();
  CBA = p["CBA"].as<quantity>();
}
} // namespace cg::nat_ang