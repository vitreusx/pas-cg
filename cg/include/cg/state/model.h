#pragma once
#include <cg/state/box.h>
#include <cg/types/amino_acid.h>
#include <cg/types/vec3.h>
#include <cg/utils/random.h>
#include <list>
#include <memory>
#include <optional>
#include <vector>

namespace cg {
class model {
public:
  model() = default;
  model(model const &other);
  model &operator=(model const &other);

  friend model operator+(model const &m1, model const &m2);
  model &operator+=(model const &m2);

  void morph_into_saw(rand_gen &gen, std::optional<double> res_bond_length,
                      double base_res_dens, bool infer_box);

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

  enum contact_type {
    UNKNOWN,
    BACK_BACK,
    BACK_SIDE,
    SIDE_BACK,
    SIDE_SIDE,
    NAT_SS
  };

  struct contact {
    residue *res1, *res2;
    double length;
    contact_type type;
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
};
} // namespace cg
