#include <cg/ckpt/parameters.h>


namespace cg::ckpt {
void parameters::load(const ioxx::xyaml::node &node) {
  node["enabled"] >> enabled;
  node["path format"] >> path_fmt;
  node["save every"] >> every;
}
} // namespace cg::ckpt