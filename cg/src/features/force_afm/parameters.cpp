#include "features/force_afm/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::fafm;

void parameters::connect(ioxx::xyaml_proxy &p) {
  enabled = p["enabled"].as<bool>();

  for (auto const &afm_p : p["AFM tips"]) {
    auto afm_proxy = p(static_cast<YAML::Node const&>(afm_p));
    auto &tip = afm_tips.emplace_back();
    tip.pulled_idx = afm_proxy["pulled idx"].as<int>();
    tip.F = afm_proxy["F"].as<vec3<quantity>>();
  }
}