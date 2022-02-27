#include "qa/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::qa;

void parameters::load(ioxx::xyaml::node const &p) {
  enabled = p["enabled"].as<bool>();
  phase_dur = p["phase duration"].as<quantity>();
  breaking_factor = p["breaking factor"].as<double>();
  min_cos_hr = p["min |cos(h, r)|"].as<double>();
  min_cos_hh = p["min |cos(h, h)| for bb"].as<double>();
  max_cos_nr = p["max cos(n, r)"].as<double>();
}