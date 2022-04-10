#include <cg/pauli/parameters.h>
#include <cg/utils/ioxx_interop.h>
namespace cg::pauli {

void parameters::link(ioxx::xyaml::proxy &proxy) {
  proxy["enabled"] & enabled;
  proxy["exclusion radius"] & r_excl;
  proxy["depth"] & depth;
}
} // namespace cg::pauli