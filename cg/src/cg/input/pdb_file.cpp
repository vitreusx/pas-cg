#include <cg/input/fields.h>
#include <cg/input/pdb_file.h>
#include <cg/input/records.h>
#include <cg/utils/quantity.h>
#include <map>
#include <set>
#include <sstream>
#include <string_view>

namespace cg {
pdb_file::pdb_file(std::istream &&is) { load(is); }

void pdb_file::load(std::istream &is) {
  std::unordered_map<int, std::unordered_map<char, bool>> ter_found;
  int cur_model_serial = 1;

  primary_model_serial = -1;

  cryst1 = decltype(cryst1)::Zero();

  for (std::string line; std::getline(is, line);) {
    if (auto opt_record = records::record::try_parse(line); opt_record) {
      auto &record_r = opt_record.value();

      if (auto model_r = record_r.cast<records::model>(); model_r) {
        cur_model_serial = model_r->serial;
      } else if (auto atom_r = record_r.cast<records::atom>(); atom_r) {
        if (atom_r->atom_name[0] == 'H')
          continue;

        auto &model = find_or_add_model(cur_model_serial);

        auto &chain = find_or_add_chain(model, atom_r->chain_id);

        auto &res =
            find_or_add_res(chain, atom_r->res_seq_num, atom_r->residue_name,
                            ter_found[cur_model_serial][chain.chain_id]);

        auto &atm = chain.atoms[atom_r->serial];

        atm.name = atom_r->atom_name;
        atm.serial = atom_r->serial;
        atm.pos = atom_r->pos * quantity("A");
        atm.parent_res = &res;

        res.atoms.push_back(&atm);

        if (!ter_found[cur_model_serial][atom_r->chain_id])
          chain.ter_serial = atom_r->serial + 1;

        if (primary_model_serial < 0)
          primary_model_serial = cur_model_serial;

      } else if (auto ssbond_r = record_r.cast<records::ssbond>(); ssbond_r) {
        disulfide_bond ss;
        ss.serial = ssbond_r->serial;
        ss.length = ssbond_r->length * quantity("A");

        auto &m = *find_model(cur_model_serial);

        auto &chain1 = *find_chain(m, ssbond_r->res[0].chain_id);
        auto &res1 = *find_res(chain1, ssbond_r->res[0].res_seq_num);
        ss.a1 = res1.find_by_name("SG");

        auto &chain2 = *find_chain(m, ssbond_r->res[1].chain_id);
        auto &res2 = *find_res(chain2, ssbond_r->res[1].res_seq_num);
        ss.a2 = res2.find_by_name("SG");

        disulfide_bonds[ss.serial] = ss;
      } else if (auto link_r = record_r.cast<records::link>(); link_r) {
        link _link;
        _link.length = link_r->length * quantity("A");

        auto &m = *find_model(cur_model_serial);

        auto &chain1 = *find_chain(m, link_r->res[0].chain_id);
        auto &res1 = *find_res(chain1, link_r->res[0].res_seq_num);
        _link.a1 = res1.find_by_name(link_r->res[0].atom_name);

        auto &chain2 = *find_chain(m, link_r->res[1].chain_id);
        auto &res2 = *find_res(chain2, link_r->res[1].res_seq_num);
        _link.a2 = res2.find_by_name(link_r->res[1].atom_name);

        links.push_back(_link);
      } else if (auto cryst1_r = record_r.cast<records::cryst1>(); cryst1_r) {
        cryst1 = cryst1_r->cell * quantity("A");
      } else if (auto ter_r = record_r.cast<records::ter>(); ter_r) {
        ter_found[cur_model_serial][ter_r->chain_id] = true;
        auto &m = *find_model(cur_model_serial);
        auto &c = *find_chain(m, ter_r->chain_id);
        c.ter_serial = ter_r->serial;
      } else if (auto end_r = record_r.cast<records::end>(); end_r) {
        return;
      }
    }
  }
}

pdb_file::pdb_file(const input::model &xmd_model) {
  char chain_id = 'A';
  size_t res_seq_num = 1, atom_serial = 1;

  std::unordered_map<input::model::residue const *, residue *> res_map;

  auto &m = find_or_add_model(primary_model_serial);

  for (auto const &xmd_chain : xmd_model.chains) {
    auto &pdb_chain = m.chains[chain_id];
    pdb_chain.chain_id = chain_id;

    for (auto const &xmd_res_ref : xmd_chain->residues) {
      auto &xmd_res = *xmd_res_ref;

      auto &ca_atom = pdb_chain.atoms[atom_serial];
      ca_atom.name = "CA";
      ca_atom.serial = atom_serial;
      ca_atom.pos = xmd_res.pos;

      auto &pdb_res = pdb_chain.residues[res_seq_num];
      pdb_res.parent_chain = &pdb_chain;
      pdb_res.seq_num = res_seq_num;
      pdb_res.name = xmd_res.type.name();

      ca_atom.parent_res = &pdb_res;
      pdb_res.atoms.push_back(&ca_atom);
      pdb_chain.order.push_back(&pdb_res);
      res_map[&xmd_res] = &pdb_res;

      ++res_seq_num;
      ++atom_serial;
    }

    ++chain_id;
  }

  size_t ss_serial = 1;
  for (auto const &xmd_cont : xmd_model.contacts) {
    auto *pdb_res1 = res_map[&*xmd_cont.res1];
    auto *pdb_res2 = res_map[&*xmd_cont.res2];

    if (xmd_cont.type == nat_cont::type::SSBOND) {
      auto &pdb_ss = disulfide_bonds[ss_serial];
      pdb_ss.a1 = pdb_res1->atoms[0];
      pdb_ss.a2 = pdb_res2->atoms[0];
      pdb_ss.serial = ss_serial;
      pdb_ss.length = xmd_cont.length;
      ++ss_serial;
    } else {
      auto &pdb_link = links.emplace_back();
      pdb_link.a1 = pdb_res1->atoms[0];
      pdb_link.a2 = pdb_res2->atoms[0];
      pdb_link.length = xmd_cont.length;
    }
  }

  cryst1 = xmd_model.model_box.cell;
}

pdb_file::atom *pdb_file::residue::find_by_name(const std::string &name) const {
  for (auto *atm : atoms) {
    if (atm->name == name)
      return atm;
  }
  return nullptr;
}

pdb_file::model *pdb_file::find_model(int model_serial) {
  auto model_iter = models.find(model_serial);
  if (model_iter != models.end())
    return &model_iter->second;
  else
    return nullptr;
}

pdb_file::model &pdb_file::find_or_add_model(int model_serial) {
  auto model_iter = models.find(model_serial);
  if (model_iter == models.end()) {
    if (models.empty())
      primary_model_serial = model_serial;

    auto &m = models[model_serial];
    m.model_serial = model_serial;

    model_iter = models.find(model_serial);
  }

  return model_iter->second;
}

pdb_file::chain *pdb_file::find_chain(model &m, char chain_id) {
  auto chain_iter = m.chains.find(chain_id);
  if (chain_iter != m.chains.end())
    return &chain_iter->second;
  else
    return nullptr;
}

pdb_file::chain &pdb_file::find_or_add_chain(model &m, char chain_id) {
  auto chain_iter = m.chains.find(chain_id);
  if (chain_iter == m.chains.end()) {
    auto &c = m.chains[chain_id];
    c.chain_id = chain_id;

    chain_iter = m.chains.find(chain_id);
  }

  return chain_iter->second;
}

pdb_file::residue *pdb_file::find_res(pdb_file::chain &c, size_t seq_num) {
  auto res_iter = c.residues.find(seq_num);
  if (res_iter != c.residues.end())
    return &res_iter->second;
  else
    return nullptr;
}

pdb_file::residue &pdb_file::find_or_add_res(chain &c, size_t seq_num,
                                             const std::string &name,
                                             bool chain_terminated) {
  auto res_iter = c.residues.find(seq_num);
  if (res_iter == c.residues.end()) {
    auto &r = c.residues[seq_num];
    r.seq_num = seq_num;
    r.name = name;
    r.parent_chain = &c;

    if (!chain_terminated)
      c.order.push_back(&r);

    res_iter = c.residues.find(seq_num);
  }

  return res_iter->second;
}

std::ostream &operator<<(std::ostream &os, const pdb_file &p) {
  bool first = true;
  for (auto const &[model_serial, model] : p.models) {
    if (!first)
      os << '\n';
    os << model;
    first = false;
  }

  if (!first)
    os << '\n';
  os << records::end().write();

  return os;
}

std::ostream &operator<<(std::ostream &os, pdb_file::model const &model) {
  records::model model_r;
  model_r.serial = model.model_serial;
  os << model_r.write();

  for (auto const &[chain_id, chain] : model.chains)
    os << '\n' << chain;

  os << '\n' << records::endmdl().write();
  return os;
}

std::ostream &operator<<(std::ostream &os, pdb_file::chain const &chain) {
  bool first = true;
  for (auto const &res : chain.order) {
    for (auto const &atm : res->atoms) {
      records::atom atom_r;
      atom_r.serial = atm->serial;
      atom_r.chain_id = chain.chain_id;
      atom_r.res_seq_num = res->seq_num;
      atom_r.residue_name = res->name;
      atom_r.atom_name = atm->name;
      atom_r.pos = atm->pos / quantity("A");

      if (!first)
        os << '\n';
      os << atom_r.write();
      first = false;
    }
  }

  auto final_res = chain.order.back();
  records::ter ter_r;
  ter_r.chain_id = chain.chain_id;
  ter_r.res_seq_num = final_res->seq_num;
  ter_r.serial = chain.ter_serial;
  ter_r.res_name = final_res->name;

  if (!first)
    os << '\n';
  os << ter_r.write();
  return os;
}

input::model pdb_file::to_model() const {
  input::model xmd_model;

  std::unordered_map<residue const *, input::model::residue *> res_map;

  auto const &m = primary_model();

  int chain_idx = 0;
  for (auto const &[chain_id, pdb_chain] : m.chains) {
    auto &xmd_chain =
        xmd_model.chains.emplace_back(std::make_unique<input::model::chain>());
    xmd_chain->chain_idx = chain_idx++;

    int res_seq_idx = 0;
    for (auto const *pdb_res : pdb_chain.order) {
      auto &xmd_res = xmd_model.residues.emplace_back(
          std::make_unique<input::model::residue>());

      xmd_res->type = amino_acid(pdb_res->name);
      xmd_res->pos = pdb_res->find_by_name("CA")->pos;
      xmd_res->parent = xmd_chain.get();
      xmd_res->seq_idx = res_seq_idx++;

      xmd_chain->residues.push_back(&*xmd_res);
      res_map[pdb_res] = &*xmd_res;
    }

    for (size_t idx = 0; idx + 1 < xmd_chain->residues.size(); ++idx) {
      auto *res1 = xmd_chain->residues[idx];
      auto *res2 = xmd_chain->residues[idx + 1];

      auto r1 = res1->pos, r2 = res2->pos;

      auto &xmd_tether = xmd_model.tethers.emplace_back();
      xmd_tether.res1 = res1;
      xmd_tether.res2 = res2;
      xmd_tether.length = norm(r2 - r1);

      if (idx + 2 < xmd_chain->residues.size()) {
        auto *res3 = xmd_chain->residues[idx + 2];
        auto r3 = res3->pos;
        auto r12_u = unit(r2 - r1), r23_u = unit(r3 - r2);
        auto theta = acos(dot(r12_u, r23_u));

        auto &xmd_angle = xmd_model.angles.emplace_back();
        xmd_angle.res1 = res1;
        xmd_angle.res2 = res2;
        xmd_angle.res3 = res3;
        xmd_angle.theta = theta;

        if (idx + 3 < xmd_chain->residues.size()) {
          auto *res4 = xmd_chain->residues[idx + 3];
          auto r4 = res4->pos;

          auto r12 = r2 - r1, r23 = r3 - r2, r34 = r4 - r3;
          auto x12_23 = cross(r12, r23), x23_34 = cross(r23, r34);
          auto x12_23_u = unit(x12_23), x23_34_u = unit(x23_34);
          auto phi = acos(dot(x12_23_u, x23_34_u));
          if (dot(x12_23, r34) < 0.0f)
            phi = -phi;

          auto &xmd_dihedral = xmd_model.dihedrals.emplace_back();
          xmd_dihedral.res1 = res1;
          xmd_dihedral.res2 = res2;
          xmd_dihedral.res3 = res3;
          xmd_dihedral.res4 = res4;
          xmd_dihedral.phi = phi;
        }
      }
    }
  }

  std::set<std::pair<residue *, residue *>> cont_res_pairs;

  for (auto const &[ss_serial, pdb_ss] : disulfide_bonds) {
    auto *res1 = pdb_ss.a1->parent_res;
    auto *res2 = pdb_ss.a2->parent_res;
    if (res1 >= res2)
      std::swap(res1, res2);
    if (cont_res_pairs.count(std::make_pair(res1, res2)) > 0)
      continue;

    auto &xmd_ss = xmd_model.contacts.emplace_back();
    xmd_ss.res1 = res_map[res1];
    xmd_ss.res2 = res_map[res2];
    xmd_ss.length = norm(xmd_ss.res1->pos - xmd_ss.res2->pos);
    xmd_ss.type = nat_cont::type::SSBOND;

    cont_res_pairs.insert(std::make_pair(res1, res2));
  }

  for (auto const &pdb_link : links) {
    auto *res1 = pdb_link.a1->parent_res;
    auto *res2 = pdb_link.a2->parent_res;
    if (res1 >= res2)
      std::swap(res1, res2);
    if (cont_res_pairs.count(std::make_pair(res1, res2)) > 0)
      continue;

    auto &xmd_cont = xmd_model.contacts.emplace_back();
    xmd_cont.res1 = res_map[res1];
    xmd_cont.res2 = res_map[res2];
    xmd_cont.length = norm(xmd_cont.res1->pos - xmd_cont.res2->pos);

    auto back1 = pdb_link.a1->in_backbone();
    auto back2 = pdb_link.a2->in_backbone();

    nat_cont::type type;
    if (back1 && back2)
      type = nat_cont::type::BACK_BACK;
    else if (back1 && !back2)
      type = nat_cont::type::BACK_SIDE;
    else if (!back1 && back2)
      type = nat_cont::type::SIDE_BACK;
    else
      type = nat_cont::type::SIDE_SIDE;

    xmd_cont.type = type;

    cont_res_pairs.insert(std::make_pair(res1, res2));
  }

  xmd_model.model_box.set_cell(cryst1);

  return xmd_model;
}

void pdb_file::add_contacts(amino_acid_data const &data, bool all_atoms) {
  std::set<std::pair<atom *, atom *>> linked_atoms;

  for (auto const &pdb_link : links) {
    linked_atoms.emplace(pdb_link.a1, pdb_link.a2);
  }

  for (auto const &[ss_serial, pdb_ss] : disulfide_bonds) {
    linked_atoms.emplace(pdb_ss.a1, pdb_ss.a2);
  }

  auto &m = primary_model();

  double alpha = pow(26.0 / 7.0, 1.0 / 6.0);
  for (auto &[chain1_id, chain1] : m.chains) {
    for (auto &[atom1_serial, atom1] : chain1.atoms) {
      if (!all_atoms && atom1.name != "CA")
        continue;

      auto res1 = amino_acid(atom1.parent_res->name);
      auto radius1 = all_atoms ? data[res1].atoms.at(atom1.name).radius
                               : data[res1].radius;
      auto seq1 = (int)atom1.parent_res->seq_num;

      for (auto &[chain2_id, chain2] : m.chains) {
        for (auto &[atom2_serial, atom2] : chain2.atoms) {
          if (!all_atoms && atom2.name != "CA")
            continue;
          if (atom1.parent_res == atom2.parent_res)
            continue;

          auto res2 = amino_acid(atom2.parent_res->name);
          auto radius2 = all_atoms ? data[res2].atoms.at(atom2.name).radius
                                   : data[res2].radius;
          auto seq2 = (int)atom2.parent_res->seq_num;

          if (abs(seq1 - seq2) < 3)
            continue;

          auto dist = norm(atom1.pos - atom2.pos);
          auto max_overlap_dist = radius1 + radius2;
          if (dist <= max_overlap_dist * alpha) {
            auto &pdb_link = links.emplace_back();
            pdb_link.a1 = &atom1;
            pdb_link.a2 = &atom2;
            pdb_link.length = dist;

            linked_atoms.emplace(&atom1, &atom2);
          }
        }
      }
    }
  }
}

pdb_file::model::model(const model &other) { *this = other; }

pdb_file::model &pdb_file::model::operator=(const model &other) {
  std::unordered_map<atom const *, atom *> atom_map;
  std::unordered_map<residue const *, residue *> res_map;

  model_serial = other.model_serial;
  chains = other.chains;

  for (auto const &[chain_idx, other_chain] : other.chains) {
    auto &cur_chain = chains[chain_idx];

    for (auto const &[atom_ser, other_atom] : other_chain.atoms) {
      auto &cur_atom = cur_chain.atoms.at(atom_ser);
      atom_map[&other_atom] = &cur_atom;
    }

    for (auto const &[seq_num, other_res] : other_chain.residues) {
      auto &cur_res = cur_chain.residues.at(seq_num);
      res_map[&other_res] = &cur_res;
    }

    for (auto &[atom_ser, atom] : cur_chain.atoms) {
      atom.parent_res = res_map.at(atom.parent_res);
    }

    for (auto &[seq_num, res] : cur_chain.residues) {
      res.parent_chain = &cur_chain;
      for (auto &atom_ptr : res.atoms) {
        atom_ptr = atom_map.at(atom_ptr);
      }
    }

    for (auto &res_ptr : cur_chain.order) {
      res_ptr = res_map.at(res_ptr);
    }
  }

  return *this;
}

pdb_file::pdb_file(const pdb_file &other) { *this = other; }

pdb_file &pdb_file::operator=(const pdb_file &other) {
  models = other.models;

  std::unordered_map<atom const *, atom *> atom_map;
  for (auto const &[model_serial, other_model] : other.models) {
    auto &cur_model = models[model_serial];
    for (auto const &[chain_id, other_chain] : other_model.chains) {
      auto &cur_chain = cur_model.chains[chain_id];
      for (auto const &[atom_serial, other_atom] : other_chain.atoms) {
        auto &cur_atom = cur_chain.atoms[atom_serial];
        atom_map[&other_atom] = &cur_atom;
      }
    }
  }

  disulfide_bonds = other.disulfide_bonds;
  for (auto &[serial, ss] : disulfide_bonds) {
    ss.a1 = atom_map[ss.a1];
    ss.a2 = atom_map[ss.a2];
  }

  links = other.links;
  for (auto &xmd_link : links) {
    xmd_link.a1 = atom_map[xmd_link.a1];
    xmd_link.a2 = atom_map[xmd_link.a2];
  }

  primary_model_serial = other.primary_model_serial;
  cryst1 = other.cryst1;

  return *this;
}

bool pdb_file::atom::in_backbone() const {
  for (auto const &back_atom : {"N", "CA", "C", "O", "OXT"}) {
    if (name == back_atom)
      return true;
  }
  return false;
}

void pdb_file::load(ioxx::xyaml::node const &node) {
  using namespace ioxx::xyaml;

  std::stringstream pdb_ss;
  pdb_ss << node.as<file>().fetch();
  load(pdb_ss);
}

pdb_file::model &pdb_file::primary_model() {
  return models.at(primary_model_serial);
}

pdb_file::model const &pdb_file::primary_model() const {
  return models.at(primary_model_serial);
}
} // namespace cg