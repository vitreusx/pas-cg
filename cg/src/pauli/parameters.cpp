#include "pauli/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::pauli;

void parameters::link(ioxx::xyaml::proxy &proxy) {
  proxy["enabled"] & enabled;
  proxy["exclusion radius"] & r_excl;
  proxy["depth"] & depth;
}