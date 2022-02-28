#pragma once
#include <ioxx/ioxx.h>

namespace cg::out {
struct report_state {
  bool first_time;
  int ord;
  std::filesystem::path output_dir;
  ioxx::xyaml::node general, current;
};
class hook {
public:
  virtual ~hook() = default;
  virtual void report_to(report_state &report) const = 0;
};
} // namespace cg::out