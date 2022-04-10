#include <cg/langevin/parameters.h>
#include <cg/utils/ioxx_interop.h>
namespace cg::lang {

void parameters::load(ioxx::xyaml::node const &p) {
  enabled = p["enabled"].as<bool>();
  gamma = p["gamma factor"].as<quantity>();
  temperature = p["temperature"].as<quantity>();
  dt = p["dt"].as<quantity>();

  if (p["type"]) {
    auto name = p["type"].as<std::string>();
    if (name == "legacy")
      type = lang_type::LEGACY;
    else if (name == "normal")
      type = lang_type::NORMAL;
  }
}
} // namespace cg::lang