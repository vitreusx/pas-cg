#include "features/pauli/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::pauli;

void parameters::connect(ioxx::xyaml_proxy &proxy) {
  proxy["enabled"] & enabled;
  proxy["exclusion radius"] & r_excl;
  proxy["depth"] & depth;
}