#include <cg/pbar/parameters.h>

namespace cg::pbar {

void parameters::link(ioxx::xyaml::proxy &proxy) {
  proxy["enabled"] & enabled;
  proxy["width"] & width;
  proxy["update period"] & update_period;
}
} // namespace cg::pbar