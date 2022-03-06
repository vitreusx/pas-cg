#pragma once
#include "morph_into_saw.h"
#include <cg/amino/amino_acid.h>
#include <cg/nat_cont/type.h>
#include <cg/types/box.h>
#include <cg/types/vec3.h>
#include <cg/utils/random.h>
#include <list>
#include <memory>
#include <optional>
#include <vector>

namespace cg::input {
class model {
public:
  model() = default;
  model(model const &other);
  model &operator=(model const &other);

  friend model operator+(model const &m1, model const &m2);
  model &operator+=(model const &m2);

  void morph_into_saw(rand_gen &gen, input::morph_into_saw_t const &params);

public:
  struct residue;
  struct chain;

  struct residue {
    chain *parent;
    int seq_idx;
    amino_acid type;
    vec3<double> pos;
  };
  std::vector<std::unique_ptr<residue>> residues;

  struct chain {
    int chain_idx;
    std::vector<residue *> residues;
  };
  std::vector<std::unique_ptr<chain>> chains;

  struct contact {
    residue *res1, *res2;
    double length;
    nat_cont::type type;
  };
  std::vector<contact> contacts;

  struct tether {
    residue *res1, *res2;
    std::optional<double> length;
  };
  std::vector<tether> tethers;

  struct angle {
    residue *res1, *res2, *res3;
    std::optional<double> theta;
  };
  std::vector<angle> angles;

  struct dihedral {
    residue *res1, *res2, *res3, *res4;
    std::optional<double> phi;
  };
  std::vector<dihedral> dihedrals;

  box<double> model_box;

  using res_map_t = std::unordered_map<residue *, int>;
};
} // namespace cg::input
