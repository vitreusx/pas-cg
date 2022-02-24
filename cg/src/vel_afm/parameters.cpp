#include "vel_afm/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::vafm;

void parameters::connect(ioxx::xyaml_proxy &p) {
  enabled = p["enabled"].as<bool>();
  H1 = p["H1"].as<quantity>();
  H2 = p["H2"].as<quantity>();

  for (auto const &afm_p : p["AFM tips"]) {
    auto afm_proxy = p(static_cast<YAML::Node const&>(afm_p));
    auto &tip = afm_tips.emplace_back();
    tip.pulled_idx = afm_proxy["pulled idx"].as<int>();
    if (auto origin_p = afm_proxy["origin"]; origin_p)
      tip.origin = origin_p.as<vec3<quantity>>();
    tip.v = afm_proxy["v"].as<vec3<quantity>>();
  }
}