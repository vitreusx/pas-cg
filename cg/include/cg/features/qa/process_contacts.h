#pragma once
#include "contact.h"
#include "free_pair.h"
#include "lj_variants.h"
#include "sync_data.h"
#include <cg/types/box.h>

namespace cg::qa {
class process_contacts {
public:
  lj_variants ljs;
  real cycle_time, cycle_time_inv, breaking_factor;
  real factor;

public:
  nitro::vector<vec3r> const *r;
  nitro::vector<vec3r> *F;
  box<real> const *box;
  nitro::set<contact> *contacts;
  real *V, *t;
  nitro::vector<sync_data> *sync;
  nitro::set<free_pair> *free_pairs;

public:
  void iter(int idx) const;
 
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::qa