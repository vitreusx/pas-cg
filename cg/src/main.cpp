#include <cg/io/seq_file.h>
#include <iostream>

int main() {
  using namespace cg;
  auto seqfile = cg::seq_file("data/glut/glut.yml");
  return 0;
}