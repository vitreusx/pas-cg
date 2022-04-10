#include <cg/pid/parameters.h>
#include <cg/utils/ioxx_interop.h>
namespace cg::pid {

void parameters::bb_t::load(ioxx::xyaml::node const &p) {
  alpha = p["alpha"].as<quantity>();
  psi_0 = p["psi_0"].as<quantity>();
  r_min = p["r_min"].as<quantity>();
  depth = p["depth"].as<quantity>();
}

void parameters::ss_t::load(ioxx::xyaml::node const &p) {
  alpha = p["alpha"].as<quantity>();
  psi_0 = p["psi_0"].as<quantity>();
}

void parameters::load(ioxx::xyaml::node const &p) {
  enabled = p["enabled"].as<bool>();
  p["include separated by 3"] >> include4;

  auto lambda_version_s = p["lambda version"].as<std::string>();
  if (lambda_version_s == "cosine")
    lambda_version = lambda_version::COSINE;
  else if (lambda_version_s == "algebraic")
    lambda_version = lambda_version::ALGEBRAIC;

  p["bb+"] >> bb_plus;
  p["bb-"] >> bb_minus;
  p["ss"] >> ss;
}
} // namespace cg::pid