#include <cg/simul/program.h>
#include <boost/archive/text_oarchive.hpp>

int main(int argc, char **argv) {
  auto prog = cg::simul::program();
  prog.main(argc, argv);
  return EXIT_SUCCESS;
}