#include <cg/nl/parameters.h>
#include <cg/utils/ioxx_interop.h>
namespace cg::nl {

void parameters::load(ioxx::xyaml::node const &node) {
  std::string algorithm_str;
  node["algorithm"] >> algorithm_str;
  if (algorithm_str == "cell")
    algorithm = CELL;
  else
    algorithm = LEGACY;
  if (node["cutoff"].as<std::string>() != "auto")
    node["cutoff"] >> cutoff;
  node["pad"] >> pad;
}
} // namespace cg::nl