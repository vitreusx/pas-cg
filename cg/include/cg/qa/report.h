#pragma once
#include "contact.h"
#include "process_contacts.h"
#include <cg/amino/sync_data.h>
#include <cg/output/hook.h>

namespace cg::qa {
class report_qa_stuff : public out::hook {
public:
  nitro::const_view<sync_data> sync_values;
  nitro::set<contact> const *contacts;
  qa::process_contacts const *process_cont;

public:
  void report_to(out::report_state &report) const override;
};
} // namespace cg::qa