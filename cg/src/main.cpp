//#include <cg/simul/program.h>
//#include <boost/archive/text_oarchive.hpp>
//
// int main(int argc, char **argv) {
//  auto prog = cg::simul::program();
//  prog.main(argc, argv);
//  return EXIT_SUCCESS;
//}

#include <iostream>
#include <ioxx/table.h>

int main() {
  using namespace ioxx::table;

  table tab = {"col1", "col2", "col3", "col4"};
  for (int row_idx = 0; row_idx < 5; ++row_idx) {
    auto &row = tab.add();
    row["col2"] = (float)row_idx;
    row["col1"] = row["col2"];
    row["col3"] = "value" + std::to_string(5 - row_idx);
    row["col4"] = false;
  }

  sl4_parser p;
  p.write(std::cout, tab);

  return EXIT_SUCCESS;
}