#include <cg/simul/default.yml.h>
#include <cg/simul/parameters.h>
namespace cg::simul {

void parameters::load(ioxx::xyaml::node const &n) {
  //  repr = YAML::Dump(n.flatten());

  n["general"] >> gen;
  n["amino acid data"] >> aa_data;
  n["input"] >> input;
  n["langevin"] >> lang;
  n["neighbor list"] >> nl;
  n["angle potentials"] >> angles;
  n["Debye-Hueckel"] >> dh;
  n["pseudo-improper dihedral"] >> pid;
  n["quasi-adiabatic"] >> qa;
  n["chirality"] >> chir;
  n["native contacts"] >> nat_cont;
  n["Pauli exclusion"] >> pauli;
  n["tether forces"] >> tether;
  n["AFM simulations"] >> afm;
  n["progress bar"] >> pbar;
  n["output"] >> out;
  n["checkpoints"] >> ckpt;
  n["local repulsive"] >> lrep;
  n["simulation box"] >> sbox;
}

ioxx::xyaml::node defaults_yml() {
  auto default_yml_str =
      std::string((const char *)default_yml, default_yml_len);
  auto default_yml_node = YAML::Load(default_yml_str.c_str());
  return ioxx::xyaml::node(default_yml_node);
}
} // namespace cg::simul