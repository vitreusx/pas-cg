#pragma once
#include "step_base.h"

namespace cg::lang {
class legacy_step : public step_base {
public:
  vect::view<vec3r> noise;

public:
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::lang