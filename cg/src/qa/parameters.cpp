#include "qa/parameters.h"
#include "utils/ioxx_interop.h"
namespace cg::qa {

void parameters::ss_spec_crit_t::load(const ioxx::xyaml::node &node) {
  node["enabled"] >> enabled;
  node["max neigh count"] >> max_neigh_count;
  node["neigh radius"] >> neigh_radius;
  node["def bond dist"] >> def_dist;
  node["max bond dist deviation"] >> max_dist_dev;
}

void parameters::disulfide_t::load(const ioxx::xyaml::node &node) {
  node["force"] >> force;
  node["special criteria"] >> spec_crit;
}

void parameters::load(ioxx::xyaml::node const &p) {
  p["enabled"] >> enabled;
  p["include separated by 3"] >> include4;
  p["phase duration"] >> phase_dur;
  p["breaking factor"] >> breaking_factor;
  p["min |cos(h, r)|"] >> min_cos_hr;
  p["min |cos(h, h)| for bb"] >> min_cos_hh;
  p["max cos(n, r)"] >> max_cos_nr;
  p["disulfides"] >> disulfide;
  p["formation tolerance"] >> formation_tolerance;
}
} // namespace cg::qa