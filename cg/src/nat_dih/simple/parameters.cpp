#include "nat_dih/simple/parameters.h"
#include "utils/ioxx_interop.h"
namespace cg::snd {

void parameters::load(ioxx::xyaml::node const &p) {
  enabled = p["enabled"].as<bool>();
  CDH = p["CDH"].as<quantity>();
}
} // namespace cg::snd