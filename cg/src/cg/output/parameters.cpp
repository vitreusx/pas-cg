#include <cg/output/parameters.h>

namespace cg::out {

void parameters::load(ioxx::xyaml::node const &from) {
  from["enabled"] >> enabled;
  if (from["emit stats every"].Scalar() != "never")
    from["emit stats every"] >> stats_every;
  if (from["emit structure every"].Scalar() != "never")
    from["emit structure every"] >> struct_every;
  from["file prefix"] >> prefix;
}
} // namespace cg::out