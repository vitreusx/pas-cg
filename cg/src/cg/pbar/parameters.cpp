#include <cg/pbar/parameters.h>
#include <cg/utils/ioxx_interop.h>
namespace cg::pbar {

void parameters::link(ioxx::xyaml::proxy &proxy) {
  proxy["enabled"] & enabled;
  proxy["width"] & width;
  proxy["update period"] & update_period;
}
} // namespace cg::pbar