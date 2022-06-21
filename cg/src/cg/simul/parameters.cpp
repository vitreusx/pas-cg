#include <cg/simul/default.yml.h>
#include <cg/simul/parameters.h>
namespace cg::simul {

void parameters::load(ioxx::xyaml::node const &n) {
  repr = YAML::Dump(n.flatten());

  n["general"] >> gen;
  n["amino acid data"] >> aa_data;
  n["input"] >> input;
  n["langevin"] >> lang;
  n["neighbor list"] >> nl;
  n["heurestic angles"] >> heur_ang;
  n["native angles"] >> nat_ang;
  n["complex native dihedrals"] >> cnd;
  n["simple native dihedrals"] >> snd;
  n["heurestic dihedrals"] >> heur_dih;
  n["constant Debye-Hueckel"] >> const_dh;
  n["relative Debye-Hueckel"] >> rel_dh;
  n["lj force variants"] >> lj_vars;
  n["pseudo-improper dihedral"] >> pid;
  n["quasi-adiabatic"] >> qa;
  n["chirality"] >> chir;
  n["native contacts"] >> nat_cont;
  n["Pauli exclusion"] >> pauli;
  n["tether forces"] >> tether;
  n["AFM"] >> afm;
  n["progress bar"] >> pbar;
  n["output"] >> out;
  n["checkpoints"] >> ckpt;
  n["local repulsive"] >> lrep;

  using prog_mode = gen::parameters::prog_mode;
  if (gen.mode == prog_mode::check_determinism) {
    out.enabled = false;
    pbar.enabled = false;
  }

  if (qa.enabled || pid.enabled)
    pauli.enabled = false;
}

ioxx::xyaml::node defaults_yml() {
  auto default_yml_str =
      std::string((const char *)default_yml, default_yml_len);
  auto default_yml_node = YAML::Load(default_yml_str.c_str());
  return ioxx::xyaml::node(default_yml_node);
}
} // namespace cg::simul