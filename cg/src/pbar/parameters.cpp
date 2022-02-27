#include "pbar/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::pbar;

void parameters::link(ioxx::xyaml::proxy &proxy) {
  proxy["enabled"] & enabled;
  proxy["width"] & width;
  proxy["update period"] & update_period;
}