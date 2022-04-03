#include <boost/serialization/access.hpp>
#include <cg/simul/program.h>

int main(int argc, char **argv) {
  auto prog = cg::simul::program();
  prog.main(argc, argv);
  return EXIT_SUCCESS;
}