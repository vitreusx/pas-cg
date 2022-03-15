#pragma once
#include "report_state.h"

namespace cg::out {
class hook {
public:
  virtual ~hook() = default;
  virtual void report_to(report_state &report) const = 0;
};
} // namespace cg::out