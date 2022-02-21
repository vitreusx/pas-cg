#include "features/pid/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::pid;

void parameters::bb_t::connect(ioxx::xyaml_proxy &p) {
  alpha = p["alpha"].as<quantity>();
  psi_0 = p["psi_0"].as<quantity>();
  r_min = p["r_min"].as<quantity>();
  depth = p["depth"].as<quantity>();
}

void parameters::ss_t::connect(ioxx::xyaml_proxy &p) {
  alpha = p["alpha"].as<quantity>();
  psi_0 = p["psi_0"].as<quantity>();
}

void parameters::connect(ioxx::xyaml_proxy &p) {
  enabled = p["enabled"].as<bool>();

  auto lambda_version_s = p["lambda version"].as<std::string>();
  if (lambda_version_s == "cosine")
    lambda_version = lambda_version_t::COSINE;
  else if (lambda_version_s == "algebraic")
    lambda_version = lambda_version_t::ALGEBRAIC;
  else
    lambda_version = lambda_version_t::UNKNOWN;

  p["bb+"] & bb_plus;
  p["bb-"] & bb_minus;
  p["ss"] & ss;
}