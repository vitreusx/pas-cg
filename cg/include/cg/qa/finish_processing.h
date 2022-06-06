#pragma once
#include "candidate.h"
#include "contact.h"
#include "free_pair.h"
#include <cg/amino/sync_data.h>

namespace cg::qa {
class finish_processing {
public:
  bool disulfide_special_criteria;

public:
  vect::vector<candidate> *candidates;
  vect::view<sync_data> sync;
  real const *t;
  vect::set<contact> *contacts;
  vect::set<free_pair> *free_pairs;
  vect::view<bool> part_of_ssbond;
  vect::vector<int> *removed;
  int *num_contacts;

public:
  void operator()() const;
};
}