#pragma once
#include "step_base.h"

namespace cg::lang {
class legacy_step : public step_base {
public:
  void operator()() const;
};
} // namespace cg::lang