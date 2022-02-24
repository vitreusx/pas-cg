#pragma once
#include "contact.h"
#include "free_pair.h"
#include "lj_variants.h"
#include <cg/amino/sync_data.h>
#include <cg/types/box.h>

namespace cg::qa {
class process_contacts {
public:
  lj_variants ljs;
  real cycle_time, cycle_time_inv;
  real factor;

  void set_factor(double breaking_factor);

public:
  nitro::const_view<vec3r> r;
  nitro::view<vec3r> F;
  box<real> const *box;
  nitro::set<contact> *contacts;
  real *V, *t;
  nitro::view<sync_data> sync;
  nitro::set<free_pair> *free_pairs;

public:
  void iter(int idx) const;

  void operator()() const;
  void omp_async() const;
};
} // namespace cg::qa