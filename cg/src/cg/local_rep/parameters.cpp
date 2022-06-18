#include <cg/local_rep/parameters.h>
#include <cg/utils/ioxx_interop.h>

namespace cg::local_rep {
void parameters::load(const ioxx::xyaml::node &n) {
  n["enabled"] >> enabled;
  n["depth"] >> depth;
}
} // namespace cg::local_rep