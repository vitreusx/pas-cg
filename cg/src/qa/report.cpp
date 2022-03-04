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

void report_qa_stuff::report_to(out::report_state &report) const {
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
  report.current["sync values"] = sync_values_file;

  ioxx::xyaml::csv<qa_contact_row> contacts_file;
  sync_values_file.path = "qa_contacts.csv";
  sync_values_file.data.header = {"i1", "i2", "type", "status"};
  for (int idx = 0; idx < contacts->size(); ++idx) {
    if (!contacts->at(idx).is_vacant()) {
      auto row = qa_contact_row(contacts->at(idx).item(), process_cont);
      contacts_file.data.rows.push_back(row);
    }
  }
  report.current["QA contacts"] = contacts_file;
}