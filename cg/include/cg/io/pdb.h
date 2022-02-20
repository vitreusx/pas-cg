#pragma once
#include <cg/state/model.h>
#include <cg/types/amino_acid.h>
#include <cg/types/vec3.h>
#include <string>
#include <vector>

namespace cg {
class pdb;
class pdb_model_emitter;
class pdb_contacts_emitter;

class pdb {
public:
  pdb(pdb const &other);
  pdb &operator=(pdb const &other);

  explicit pdb(std::istream &&is);
  explicit pdb(cg::model const &xmd_model);

  friend std::ostream &operator<<(std::ostream &os, pdb const &p);
  pdb_model_emitter emit_model(int model_serial) const;
  pdb_contacts_emitter emit_contacts() const;

  void add_contacts(amino_acid_data const &data, bool all_atoms = true);

  enum class contact_deriv { NONE, FROM_ATOMS, FROM_RESIDUES };
  cg::model to_model() const;

  struct atom;
  struct residue;
  struct chain;
  struct model;

  struct atom {
    std::string name;
    size_t serial;
    vec3<double> pos;
    residue *parent_res;

    bool in_backbone() const;
  };

  struct residue {
    chain *parent_chain;
    size_t seq_num;
    std::string name;
    std::vector<atom *> atoms;

    atom *find_by_name(std::string const &name) const;
  };

  struct chain {
    char chain_id;
    std::unordered_map<size_t, atom> atoms;
    std::unordered_map<size_t, residue> residues;
    std::vector<residue *> order;
    size_t ter_serial;
  };
  std::unordered_map<char, chain> chains;

  struct disulfide_bond {
    size_t serial;
    atom *a1, *a2;
    double length;
  };
  std::unordered_map<size_t, disulfide_bond> disulfide_bonds;

  struct link {
    atom *a1, *a2;
    double length;
  };
  std::vector<link> links;

  vec3<double> cryst1;

private:
  friend class pdb_model_emitter;
  friend class pdb_contacts_emitter;

  chain *find_chain(char chain_id);
  residue *find_res(chain &c, size_t seq_num);
};

class pdb_model_emitter {
public:
  friend std::ostream &operator<<(std::ostream &os,
                                  pdb_model_emitter const &emitter);

private:
  friend class pdb;

  pdb const &owner;
  int model_serial;
  explicit pdb_model_emitter(pdb const &owner, int model_serial);
};

class pdb_contacts_emitter {
public:
  friend std::ostream &operator<<(std::ostream &os,
                                  pdb_contacts_emitter const &emitter);

private:
  friend class pdb;

  pdb const &owner;
  explicit pdb_contacts_emitter(pdb const &owner);
};
} // namespace cg