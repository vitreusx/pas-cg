#pragma once
#include "candidate.h"
#include "contact.h"
#include "free_pair.h"
#include "sync_data.h"

namespace cg::qa {
class process_candidates {
public:
  nitro::vector<candidate> *candidates;
  nitro::vector<sync_data> *sync;
  real const *t;
  nitro::set<contact> *contacts;
  nitro::set<free_pair> *free_pairs;

public:
  template <typename E> void iter(candidate_expr<E> const &candidate) const;
  void operator()() const;
};
}; // namespace cg::qa