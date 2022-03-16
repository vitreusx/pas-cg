#include "nat_cont/report.h"
#include "utils/quantity.h"
#include <ioxx/ioxx.h>
namespace cg::nat_cont {

struct nat_cont_row {
  nat_cont cont;
  report_stuff const *master;
  nat_cont_row(nat_cont const &cont, report_stuff const &master)
      : cont(cont), master{&master} {};

  void connect(ioxx::row_proxy &row) {
    row["i1"] << cont.i1();
    row["chain1"] << master->chain_idx[cont.i1()];
    row["i2"] << cont.i2();
    row["chain2"] << master->chain_idx[cont.i2()];
    row["nat_dist[A]"] << cg::quantity(cont.nat_dist()).in("A");
    row["cur_dist[A]"] << master->cur_dist(cont).in("A");
    row["is active"] << master->is_active(cont);
    row["formed once"] << cont.formed();
    if (cont.formed())
      row["formation time"] << cont.formation_t();
  }
};

enum class count_type : int { INTRA, INTER };

struct contact_count {
  int count[2] = {0, 0};

  int &operator[](count_type type) { return count[(int)type]; }

  int all() const {
    return count[(int)count_type::INTRA] + count[(int)count_type::INTER];
  }
};

cg::quantity report_stuff::cur_dist(nat_cont cont) const {
  return norm(r[cont.i2()] - r[cont.i1()]);
}

bool report_stuff::is_active(nat_cont cont) const {
  return cur_dist(cont) < cont.nat_dist() * params->active_thr;
}

void report_stuff::report_to(out::report_data &report) const {
  auto &nc_node = report.for_snap["native contacts"];

  if (report.report_files) {
    ioxx::xyaml::csv<nat_cont_row> nat_conts_file;
    nat_conts_file.path = "nat_conts.csv";
    nat_conts_file.data.header = {"i1",        "chain1",      "i2",
                                  "chain2",    "nat_dist[A]", "cur_dist[A]",
                                  "is active", "formed once", "formation time"};

    for (int idx = 0; idx < all_contacts.size(); ++idx) {
      auto cont = all_contacts[idx];
      nat_cont_row row(cont, *this);
      nat_conts_file.data.rows.push_back(row);
    }

    nc_node["contacts"] = nat_conts_file;
  }

  contact_count num_all, num_bb, num_bs, num_ss, num_ssbond;

  for (int idx = 0; idx < all_contacts.size(); ++idx) {
    auto cont = all_contacts[idx];
    if (is_active(cont)) {
      auto count_type = chain_idx[cont.i1()] == chain_idx[cont.i2()]
                            ? count_type::INTRA
                            : count_type::INTER;

      ++num_all[count_type];
      switch (cont.type()) {
      case type::BACK_BACK:
      case type::ANY:
        ++num_bb[count_type];
        break;
      case type::BACK_SIDE:
      case type::SIDE_BACK:
        ++num_bs[count_type];
        break;
      case type::SIDE_SIDE:
        ++num_ss[count_type];
        break;
      case type::SSBOND:
        ++num_ssbond[count_type];
        break;
      }
    }
  }

  auto &active_node = nc_node["num of active contacts"];
  active_node["all"] = num_all.all();
  active_node["back-back"] = num_bb.all();
  active_node["back-back (intra)"] = num_bb[count_type::INTRA];
  active_node["back-back (inter)"] = num_bb[count_type::INTER];
  active_node["back-side"] = num_bs.all();
  active_node["back-side (intra)"] = num_bs[count_type::INTRA];
  active_node["back-side (inter)"] = num_bs[count_type::INTER];
  active_node["side-side"] = num_ss.all();
  active_node["side-side (intra)"] = num_ss[count_type::INTRA];
  active_node["side-side (inter)"] = num_ss[count_type::INTER];
}
} // namespace cg::nat_cont