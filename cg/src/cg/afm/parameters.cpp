#include <cg/afm/parameters.h>

namespace cg::afm {

void parameters::fafm_t::load(const ioxx::xyaml::node &n) {
  n["force"] >> force;
}

void parameters::vafm_t::load(const ioxx::xyaml::node &n) {
  n["velocity"] >> vel;
  n["H1"] >> H1;
  n["H2"] >> H2;
}

void parameters::pull_rel_t::load(const ioxx::xyaml::node &n) {
  n["time"] >> time;
}

void parameters::tip_params_t::load(const ioxx::xyaml::node &n) {
  n["type"] >> type;
  n["force params"] >> force_afm;
  n["velocity params"] >> vel_afm;
}

void parameters::load(ioxx::xyaml::node const &n) {
  n["perform"] >> perform;
  n["type"] >> type;
  n["pull-release params"] >> pull_rel;
  n["tip params"] >> tip_params;
}
} // namespace cg::afm