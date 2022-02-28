#include "nat_cont/report.h"
#include "utils/quantity.h"
#include <ioxx/ioxx.h>
using namespace cg::nat_cont;

struct nat_cont_row {
  nat_cont cont;
  report_stuff const *master;
  nat_cont_row(nat_cont const &cont, report_stuff const &master)
      : cont(cont), master{&master} {};

  void connect(ioxx::row_proxy &row) {
    row["i1"] << cont.i1();
    row["i2"] << cont.i2();
    row["nat_dist[A]"] << cg::quantity(cont.nat_dist()).in("A");
    auto cur_dist = norm(master->r[cont.i2()] - master->r[cont.i1()]);
    row["cur_dist[A]"] << cg::quantity(cur_dist).in("A");
    row["is active"] << (cur_dist <
                         master->params->active_thr * cont.nat_dist());
  }
};

void report_stuff::report_to(out::report_state &report) const {
  ioxx::xyaml::csv<nat_cont_row> nat_conts_file;
  nat_conts_file.path = "nat_conts.csv";
  nat_conts_file.data.header = {"i1", "i2", "nat_dist[A]", "cur_dist[A]",
                                "is active"};
  for (int idx = 0; idx < all_nat_conts.size(); ++idx) {
    nat_cont_row row(all_nat_conts[idx], *this);
    nat_conts_file.data.rows.push_back(row);
  }
  report.current["native contacts"] = nat_conts_file;
}