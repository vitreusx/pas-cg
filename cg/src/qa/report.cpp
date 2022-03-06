#include "qa/report.h"
#include "utils/quantity.h"
#include <ioxx/ioxx.h>
using namespace cg::qa;

struct sync_values_row {
  int idx;
  cg::sync_data sync;

  void connect(ioxx::row_proxy &row) {
    row["idx"] << idx;
    row["back"] << (int)sync.back();
    row["side (all)"] << (int)sync.side_all();
    row["side (hydrophobic)"] << (int)sync.side_hydrophobic();
    row["side (polar)"] << (int)sync.side_polar();
  }
};

struct qa_contact_row {
  contact cont;
  process_contacts const *proc_cont;
  explicit qa_contact_row(contact const &cont,
                          process_contacts const *proc_cont)
      : cont{cont}, proc_cont{proc_cont} {};

  static std::string contact_type_str(contact_type const &type) {
    if ((short)type == (short)contact_type::BACK_BACK())
      return "bb";
    else if ((short)type == (short)contact_type::BACK_SIDE())
      return "bs";
    else if ((short)type == (short)contact_type::SIDE_BACK())
      return "sb";
    else
      return "ss";
  }

  static std::string contact_status_str(contact_status const &status) {
    if (status == cg::qa::FORMING_OR_FORMED)
      return "forming or formed";
    else
      return "breaking";
  }

  void connect(ioxx::row_proxy &row) {
    row["i1"] << cont.i1();
    row["i2"] << cont.i2();
    row["type"] << contact_type_str(cont.type());
    row["status"] << contact_status_str(cont.status());
    row["ref_time[tau]"] << cg::quantity(cont.ref_time()).in("tau");
    row["saturation"] << proc_cont->saturation_value(cont);
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

void report_qa_stuff::report_to(out::report_state &report) const {
  auto &qa_node = report.current["quasi-adiabatic"];

  ioxx::xyaml::csv<sync_values_row> sync_values_file;
  sync_values_file.path = "sync_values.csv";
  sync_values_file.data.header = {"idx", "back", "side (all)",
                                  "side (hydrophobic)", "side (polar)"};
  for (int idx = 0; idx < sync_values.size(); ++idx) {
    sync_values_row row;
    row.idx = idx;
    row.sync = sync_values[idx];
    sync_values_file.data.rows.push_back(row);
  }
  qa_node["sync values"] = sync_values_file;

  ioxx::xyaml::csv<qa_contact_row> contacts_file;
  contacts_file.path = "qa_contacts.csv";
  contacts_file.data.header = {
      "i1", "i2", "type", "status", "ref_time[tau]", "saturation"};
  for (int idx = 0; idx < contacts->size(); ++idx) {
    if (contacts->at(idx).has_item()) {
      auto row = qa_contact_row(contacts->at(idx).item(), process_cont);
      contacts_file.data.rows.push_back(row);
    }
  }
  qa_node["contacts"] = contacts_file;

  contact_count num_all, num_bb, num_bs, num_ss, num_dyn_ss;
  for (int idx = 0; idx < contacts->size(); ++idx) {
    if (contacts->at(idx).is_vacant())
      continue;

    auto cont = contacts->at(idx).item();
    auto count_type = chain_idx[cont.i1()] == chain_idx[cont.i2()]
                          ? count_type::INTRA
                          : count_type::INTER;

    ++num_all[count_type];
    if (cont.type() == contact_type::BACK_BACK()) {
      ++num_bb[count_type];
    } else if (cont.type() == contact_type::BACK_SIDE() ||
               cont.type() == contact_type::SIDE_BACK()) {
      ++num_bs[count_type];
    } else {
      static auto ss_type = contact_type::SIDE_SIDE(aa_code::CYS, aa_code::CYS);
      if (dyn_ss && cont.type() == ss_type)
        ++num_dyn_ss[count_type];
      else
        ++num_ss[count_type];
    }
  }

  auto &active_node = qa_node["num of active contacts"];
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
  active_node["dynamic ssbonds"] = num_dyn_ss.all();
  active_node["dynamic ssbonds (intra)"] = num_dyn_ss[count_type::INTRA];
  active_node["dynamic ssbonds (inter)"] = num_dyn_ss[count_type::INTER];
}