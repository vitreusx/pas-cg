#include "dh/const/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::const_dh;

void parameters::link(ioxx::xyaml::proxy &p) {
  p["enabled"] & enabled;
  p["screening distance"] & screening_dist;
  p["permittivity"] & permittivity;
}