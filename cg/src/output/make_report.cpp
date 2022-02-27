#include "output/make_report.h"
using namespace cg::output;

void make_report::operator()() const {
  ioxx::xyaml::node report;
  report.loc = output_dir / std::to_string(*ord) / "output.yml";
  for (auto const *hook : hooks)
    hook->report_to(report);
  report.save();
}