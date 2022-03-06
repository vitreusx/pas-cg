#include "afm/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::afm;

void parameters::load(ioxx::xyaml::node const &node) {
  node["H1"] >> H1;
  node["H2"] >> H2;
  if (node["tips"]) {
    for (auto entry : node["tips"]) {
      auto xentry = node.child(entry);
      auto type_str = xentry["type"].as<std::string>();
      tip_type type = type_str == "const velocity" ? tip_type::CONST_VEL
                                                   : tip_type::CONST_FORCE;

      if (auto sin_res_node = xentry["single residue"]; sin_res_node) {
        single_res_t tip;
        tip.type = type;
        xentry["direction"] >> tip.dir;
        tip.res_idx = sin_res_node.as<int>();
        tips.push_back(tip);
      } else if (auto pulled_node = xentry["pulled-apart chain"]; pulled_node) {
        pulled_apart_t tip;
        tip.type = type;
        xentry["magnitude"] >> tip.mag;
        tip.chain_idx = pulled_node.as<int>();
        tips.push_back(tip);
      }
    }
  }
  enabled = !tips.empty();
}