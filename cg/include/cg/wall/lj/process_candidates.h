#pragma once
#include "data.h"

namespace cg::wall::lj {
class process_candidates {
public:
  vect::const_view<vec3r> r;
  vect::view<wall> walls;
  vect::set<connection> *conns;
  vect::vector<int> *removed;
  vect::vector<candidate> *candidates;
  vect::view<bool> is_connected;

public:
  void operator()() const;
};
} // namespace cg::wall::lj