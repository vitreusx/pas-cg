#pragma once
#include <ioxx/ioxx.h>

namespace cg::output {
class hook {
public:
  virtual ~hook() = default;
  virtual void report_to(ioxx::xyaml::node const &node) const = 0;
};
} // namespace cg::output