#pragma once
#include <cg/amino/aa_data.h>
#include <cg/input/model.h>
#include <cg/types/vec3.h>
#include <string>
#include <vector>

namespace cg {
class pdb_file;

struct pdb_load_options {
  bool skip_unknown = true;
  std::unordered_map<std::string, amino_acid> aliases;
  void load(ioxx::xyaml::node const &node);
};

class pdb_file {
public:
  pdb_file() = default;
  explicit pdb_file(input::model const &xmd_model);

  pdb_file(pdb_file const &other);
  pdb_file &operator=(pdb_file const &other);

  friend std::ostream &operator<<(std::ostream &os, pdb_file const &p);

  void add_contacts(amino_acid_data const &data, bool all_atoms = true);

  enum class contact_deriv {
    FROM_ATOMS,
    FROM_RESIDUES
  };
  input::model to_model() const;

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
    amino_acid type;
    std::vector<atom *> atoms;

    atom *find_by_name(std::string const &name) const;
  };

  struct chain {
    char chain_id;
    std::unordered_map<size_t, atom> atoms;
    std::unordered_map<size_t, residue> residues;
    std::vector<residue *> order;
    size_t ter_serial;

    friend std::ostream &operator<<(std::ostream &os, chain const &chain);
  };

  struct model {
    int model_serial;
    std::unordered_map<char, chain> chains;

    model() = default;
    model(model const &other);
    model &operator=(model const &other);

    friend std::ostream &operator<<(std::ostream &os, model const &model);
  };

  int primary_model_serial;
  std::map<int, model> models;
  model &primary_model();
  model const &primary_model() const;

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

  void load(ioxx::xyaml::node const &node);
  void load(std::istream &source, pdb_load_options const &load_opts);

  model *find_model(int model_serial);
  model &find_or_add_model(int model_serial);

  chain *find_chain(model &m, char chain_id);
  chain &find_or_add_chain(model &m, char chain_id);

  residue *find_res(chain &c, size_t seq_num);
  residue &find_or_add_res(chain &c, size_t seq_num, const std::string &name,
                           bool chain_terminated,
                           pdb_load_options const &load_opts);
};
} // namespace cg