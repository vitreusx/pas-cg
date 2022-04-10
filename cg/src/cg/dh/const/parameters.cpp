#include <cg/dh/const/parameters.h>
#include <cg/utils/ioxx_interop.h>
namespace cg::const_dh {

void parameters::link(ioxx::xyaml::proxy &p) {
  p["enabled"] & enabled;
  p["screening distance"] & screening_dist;
  p["permittivity"] & permittivity;
}
} // namespace cg::const_dh