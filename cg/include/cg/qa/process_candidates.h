#pragma once
#include "candidate.h"
#include "contact.h"
#include "free_pair.h"
#include <cg/amino/sync_data.h>

namespace cg::qa {
class process_candidates {
public:
  bool disulfide_special_criteria;

public:
  nitro::vector<candidate> *candidates;
  nitro::view<sync_data> sync;
  real const *t;
  nitro::set<contact> *contacts;
  nitro::set<free_pair> *free_pairs;
  nitro::view<bool> part_of_ssbond;

public:
  template <typename E> void iter(candidate_expr<E> const &candidate) const;
  void operator()() const;
};
}; // namespace cg::qa