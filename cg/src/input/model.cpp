#include "input/model.h"
#include "types/mat3x3.h"
#include <Eigen/Geometry>
#include <unordered_map>

namespace cg::input {
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

template <typename U>
static inline Eigen::Vector3<U> perp_v(Eigen::Vector3<U> const &v) {
  auto perp1 = v.cross(Eigen::Vector3<U>::UnitX());
  auto perp2 = v.cross(Eigen::Vector3<U>::UnitY());
  auto &perp = perp1.norm() > perp2.norm() ? perp1 : perp2;
  return perp.normalized();
}

void model::morph_into_saw(rand_gen &gen,
                           input::morph_into_saw_t const &params) {

  using U = double;
  auto max_spread = M_PI / 3;
  auto max_around = M_2_PI;

  auto vol = (double)residues.size() / params.residue_density;
  auto cell_a = cbrt(vol);

  double def_bond_distance = params.bond_distance.value_or(0);

  bool found_conformation = false;
  for (int retry_idx = 0; retry_idx < params.num_of_retries; ++retry_idx) {
    for (auto const &xmd_chain : chains) {
      Eigen::Vector3<U> pos{gen.uniform<U>(-cell_a / 2.0, cell_a / 2.0),
                            gen.uniform<U>(-cell_a / 2.0, cell_a / 2.0),
                            gen.uniform<U>(-cell_a / 2.0, cell_a / 2.0)};

      Eigen::Vector3<U> dir = convert<U>(gen.sphere<U>());

      for (size_t res_idx = 0; res_idx < xmd_chain->residues.size();
           ++res_idx) {
        auto next = pos;
        if (res_idx + 1 < xmd_chain->residues.size()) {
          U bond_length;
          if (params.bond_distance.has_value()) {
            bond_length = def_bond_distance;
          } else {
            auto pos1 = xmd_chain->residues[res_idx]->pos;
            auto pos2 = xmd_chain->residues[res_idx + 1]->pos;
            bond_length = norm(pos2 - pos1);
          }

          next += dir * bond_length;

          U spread_angle = gen.uniform(-max_spread, max_spread);
          auto spread = Eigen::AngleAxisd(spread_angle, perp_v(dir));

          U around_angle = gen.uniform(-max_around, max_around);
          auto around = Eigen::AngleAxisd(around_angle, dir);

          dir = (around * spread * dir).normalized();
        }

        xmd_chain->residues[res_idx]->pos = pos;
        pos = next;
      }
    }

    if (params.intersection_at != 0) {
      double ix_dist_inv = 1.0 / params.intersection_at;
      bool do_intersect = false;
      for (int idx1 = 0; idx1 < (int)residues.size(); ++idx1) {
        auto &res1 = residues[idx1];
        for (int idx2 = idx1 + 1; idx2 < (int)residues.size(); ++idx2) {
          auto &res2 = residues[idx2];
          if (norm_inv(res2->pos - res1->pos) > ix_dist_inv) {
            do_intersect = true;
            break;
          }
        }

        if (do_intersect)
          break;
      }

      if (!do_intersect) {
        found_conformation = true;
        break;
      }
    }
    else {
      found_conformation = true;
      break;
    }
  }

  if (!found_conformation)
    throw std::runtime_error("conformation not found!");

  if (params.infer_box) {
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

void model::morph_into_line(double bond_dist) {
  vec3d cur = vec3d::Zero();
  for (auto const &chain : chains) {
    for (auto *res : chain->residues) {
      res->pos = cur;
      cur.z() += bond_dist;
    }
  }
}
} // namespace cg::input
