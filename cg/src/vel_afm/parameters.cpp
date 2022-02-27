#include "vel_afm/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::vafm;

void parameters::load(ioxx::xyaml::node const &p) {
  enabled = p["enabled"].as<bool>();
  H1 = p["H1"].as<quantity>();
  H2 = p["H2"].as<quantity>();

  for (auto const &afm_p : p["AFM tips"]) {
    auto afm_tip_node = p.child(afm_p);
    auto &tip = afm_tips.emplace_back();
    afm_tip_node["pulled idx"] >> tip.pulled_idx;
    afm_tip_node["origin"] >> tip.origin;
    afm_tip_node["v"] >> tip.v;
  }
}