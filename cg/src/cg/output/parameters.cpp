#include <cg/output/parameters.h>
#include <cg/utils/ioxx_interop.h>
namespace cg::out {

void parameters::load(ioxx::xyaml::node const &from) {
  from["enabled"] >> enabled;
  from["emit stats every"] >> stats_every;
  from["emit structure every"] >> struct_every;
  from["file prefix"] >> prefix;
}
} // namespace cg::out