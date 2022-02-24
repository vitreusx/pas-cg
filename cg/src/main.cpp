#include <cg/simul/simulation.h>

int main(int argc, char **argv) {
  //  auto params_yml =
  //  ioxx::xyaml_node::from_path("data/default/inputfile.yml"); auto proxy =
  //  ioxx::xyaml_proxy(params_yml); std::cout << proxy["general"]["total
  //  time"].as<std::string>(); return EXIT_SUCCESS;
  return cg::simul::simulation()(argc, argv);
}