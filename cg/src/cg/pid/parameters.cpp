#include <cg/pid/parameters.h>

namespace cg::pid {
void parameters::lambda_params::load(const ioxx::xyaml::node &node) {
  node["alpha"] >> alpha;
  node["psi_0"] >> psi_0;
}

void parameters::load(const ioxx::xyaml::node &node) {
  node["enabled"] >> enabled;
  node["include (i, i+4)"] >> include4;

  auto lambda_n = node["lambda"];
  auto lambda_variant_name = lambda_n["variant"].as<std::string>();
  if (lambda_variant_name == "cosine")
    lambda_variant = lambda_version::COSINE;
  else if (lambda_variant_name == "algebraic")
    lambda_variant = lambda_version::ALGEBRAIC;

  lambda_n["bb+"] >> bb_plus_lambda;
  lambda_n["bb-"] >> bb_minus_lambda;
  lambda_n["ss"] >> ss_lambda;

  auto forces_n = node["forces"];
  forces_n["variant"] >> def_force_variant;
  bb_plus_force.variant = def_force_variant;
  bb_minus_force.variant = def_force_variant;
  ss_force.variant = def_force_variant;

  forces_n["bb-"] >> bb_minus_force;
  forces_n["bb+"] >> bb_plus_force;
  forces_n["ss"] >> ss_force;
}
} // namespace cg::pid