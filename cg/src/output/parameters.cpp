#include "output/parameters.h"
#include "utils/ioxx_interop.h"
namespace cg::out {

void parameters::load(ioxx::xyaml::node const &from) {
  from["enabled"] >> enabled;
  from["stats period"] >> stats_period;
  from["file period"] >> file_period;
  from["output dir"] >> output_dir;
}
} // namespace cg::out