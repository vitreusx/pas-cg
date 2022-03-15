#pragma once
#include "report_data.h"

namespace cg::out {
class hook {
public:
  virtual ~hook() = default;
  virtual void report_to(report_data &report) const = 0;
};
} // namespace cg::out