#pragma once
#include <array>
#include <cg/amino/amino_acid.h>
#include <cg/nl/data.h>
#include <cg/sbox/pbc.h>
#include <cg/vect/vect.h>

namespace cg::qa {
class legacy_impl {
public:
  real cutoff;
  real charge[amino_acid::NUM_TYPES];
  bool include4;

public:
  int jq, jq2;

  using qa_data_t = std::array<int, 4>;
  vect::view<vect::vector<qa_data_t>> kqist;

  using es_data_t = std::array<int, 3>;
  vect::vector<es_data_t> *kcist;

  vect::const_view<vec3r> r;
  sbox::pbc<real> const *simul_box;
  nl::data *nl;
  vect::const_view<int> chain_idx, seq_idx;
  vect::const_view<amino_acid> atype;

public:
  void update_verlet_list() const;
  void evalcpot() const;
};
} // namespace cg::qa