#include "nat_ang/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::nat_ang;

void parameters::load(ioxx::xyaml::node const &p) {
  enabled = p["enabled"].as<bool>();
  k = p["k"].as<quantity>();
}