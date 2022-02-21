#include <cg/features/parameters.h>

int main() {
  using namespace cg;
  auto paramfile = parameters("data/default/inputfile.yml");
  return 0;
}