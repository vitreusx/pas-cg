#include "state/model.h"
#include "types/mat3x3.h"
#include <unordered_map>

namespace cg {
model::model(const model &other) { *this = other; }

model &model::operator=(const model &other) {
  std::unordered_map<residue const *, residue *> res_map;

  residues.clear();
  for (auto const &other_res : other.residues) {
    auto &res = residues.emplace_back(std::make_unique<residue>(*other_res));

    res_map[&*other_res] = &*res;
  }

  chains.clear();
  for (auto const &other_chain : other.chains) {
    auto &xmd_chain =
        chains.emplace_back(std::make_unique<chain>(*other_chain));

    for (residue *&res : xmd_chain->residues) {
      res = res_map[res];
      res->parent = xmd_chain.get();
    }
  }

  tethers = other.tethers;
  for (auto &tether_ : tethers) {
    tether_.res1 = res_map[tether_.res1];
    tether_.res2 = res_map[tether_.res2];
  }

  contacts = other.contacts;
  for (auto &cont : contacts) {
    cont.res1 = res_map[cont.res1];
    cont.res2 = res_map[cont.res2];
  }

  angles = other.angles;
  for (auto &_angle : angles) {
    _angle.res1 = res_map[_angle.res1];
    _angle.res2 = res_map[_angle.res2];
    _angle.res3 = res_map[_angle.res3];
  }

  dihedrals = other.dihedrals;
  for (auto &_dihedral : dihedrals) {
    _dihedral.res1 = res_map[_dihedral.res1];
    _dihedral.res2 = res_map[_dihedral.res2];
    _dihedral.res3 = res_map[_dihedral.res3];
    _dihedral.res4 = res_map[_dihedral.res4];
  }

  model_box = other.model_box;

  return *this;
}

model &model::operator+=(const model &m2) {
  auto m2_copy = m2;

  for (auto &res : m2_copy.residues)
    residues.emplace_back(std::move(res));
  for (auto &chain : m2_copy.chains)
    chains.emplace_back(std::move(chain));

  tethers.insert(tethers.end(), m2_copy.tethers.begin(), m2_copy.tethers.end());
  contacts.insert(contacts.end(), m2_copy.contacts.begin(),
                  m2_copy.contacts.end());
  angles.insert(angles.end(), m2_copy.angles.begin(), m2_copy.angles.end());
  dihedrals.insert(dihedrals.end(), m2_copy.dihedrals.begin(),
                   m2_copy.dihedrals.end());

  return *this;
}

model operator+(const model &m1, const model &m2) {
  model sum;
  sum += m1;
  sum += m2;
  return sum;
}

template <typename U> static inline vec3<U> perp_v(vec3<U> const &v) {
  vec3<U> perp1 = cross(v, vec3<U>::UnitX());
  vec3<U> perp2 = cross(v, vec3<U>::UnitY());
  auto &perp = norm(perp1) > norm(perp2) ? perp1 : perp2;
  return unit(perp);
}

void model::morph_into_saw(rand_gen &gen, std::optional<double> res_bond_length,
                           double base_res_dens, bool infer_box) {

  using U = double;
  auto max_spread = M_PI / 12.0;
  auto max_around = M_PI;

  auto vol = (double)residues.size() / base_res_dens;
  auto cell_a = cbrt(vol);

  for (auto const &xmd_chain : chains) {
    vec3<U> pos{gen.uniform<U>(-cell_a / 2.0, cell_a / 2.0),
                gen.uniform<U>(-cell_a / 2.0, cell_a / 2.0),
                gen.uniform<U>(-cell_a / 2.0, cell_a / 2.0)

    };

    vec3<U> dir = gen.sphere<U>();

    for (size_t res_idx = 0; res_idx < xmd_chain->residues.size(); ++res_idx) {
      auto next = pos;
      if (res_idx + 1 < xmd_chain->residues.size()) {
        U bond_length;
        if (res_bond_length) {
          bond_length = res_bond_length.value();
        } else {
          auto pos1 = xmd_chain->residues[res_idx]->pos;
          auto pos2 = xmd_chain->residues[res_idx + 1]->pos;
          bond_length = norm(pos2 - pos1);
        }

        next += dir * bond_length;

        U spread_angle = gen.uniform(-max_spread, max_spread);
        auto spread = rot_around(perp_v(dir), spread_angle);

        U around_angle = gen.uniform(-max_around, max_around);
        auto around = rot_around(dir, around_angle);

        dir = unit(around * spread * dir);
      }

      xmd_chain->residues[res_idx]->pos = pos;
      pos = next;
    }
  }

  if (infer_box) {
    auto min_val = std::numeric_limits<double>::min();
    auto x_max = min_val, y_max = min_val, z_max = min_val;
    auto max_val = std::numeric_limits<double>::max();
    auto x_min = max_val, y_min = max_val, z_min = max_val;

    for (auto const &res : residues) {
      auto const &p = res->pos;
      x_min = std::min(x_min, p.x());
      x_max = std::max(x_max, p.x());
      y_min = std::min(y_min, p.x());
      y_max = std::max(y_max, p.x());
      z_min = std::min(z_min, p.x());
      z_max = std::max(z_max, p.x());
    }

    auto x_ext = x_max - x_min, y_ext = y_max - y_min, z_ext = z_max - z_min;

    auto x_ext_inv = 1.0 / x_ext, y_ext_inv = 1.0 / y_ext,
         z_ext_inv = 1.0 / z_ext;

    model_box.cell = vec3<U>{x_ext, y_ext, z_ext};
    model_box.cell_inv = vec3<U>{x_ext_inv, y_ext_inv, z_ext_inv};
  }
}
} // namespace cg
