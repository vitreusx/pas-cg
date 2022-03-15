#include "output/parameters.h"
#include "utils/ioxx_interop.h"
namespace cg::out {

void parameters::load(ioxx::xyaml::node const &from) {
  from["enabled"] >> enabled;
  from["period"] >> period;
  from["output dir"] >> output_dir;
}
} // namespace cg::out