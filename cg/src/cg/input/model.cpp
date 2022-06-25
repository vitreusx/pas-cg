#include <Eigen/Geometry>
#include <cg/input/model.h>
#include <cg/types/mat3x3.h>
#include <unordered_map>

namespace cg::input {
model::model(const model &other) {
  *this = other;
}

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
                           const input::morph_into_saw_t &params) {
#ifndef COMPAT_MODE
  morph_into_saw_def(gen, params);
#else
  morph_into_saw_f77(gen, params);
#endif
}

void model::morph_into_saw_def(rand_gen &gen,
                               input::morph_into_saw_t const &params) {

  using U = double;
  auto max_spread = M_PI / 3;
  auto max_around = M_2_PI;

  Eigen::AlignedBox3d box;
  if (params.init_box_density > 0) {
    auto vol = (double)residues.size() / params.init_box_density;
    auto cell_a = cbrt(vol);
    box.min() = {-cell_a / 2.0, -cell_a / 2.0, -cell_a / 2.0};
    box.max() = -box.min();
    if (params.with_pbc) {
      Eigen::Vector3d box_dim = box.max() - box.min();
      model_box.set_cell(box_dim);
    }
  } else {
    if (params.with_pbc)
      throw std::runtime_error("with_pbc requires sample_from_box");

    box.min() = Eigen::Vector3d::Zero();
    box.max() = -box.min();
  }

  double def_bond_distance = params.bond_distance.value_or(0);

  bool found_conformation = false;
  for (int retry_idx = 0; retry_idx < params.num_of_retries; ++retry_idx) {
    for (auto const &xmd_chain : chains) {
      Eigen::Vector3<U> pos{gen.uniform<U>(box.min().x(), box.max().x()),
                            gen.uniform<U>(box.min().y(), box.max().y()),
                            gen.uniform<U>(box.min().z(), box.max().z())};

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
            if (params.with_pbc)
              bond_length = norm(model_box.wrap(pos1, pos2));
            else
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
    } else {
      found_conformation = true;
      break;
    }
  }

  if (!found_conformation)
    throw std::runtime_error("conformation not found!");
}

static Eigen::Vector3d sample_from_box(Eigen::AlignedBox3d const &box,
                                       rand_gen &gen) {
  Eigen::Vector3d v;
  for (int i = 0; i < 3; ++i) {
    v[i] = (box.max()[i] - box.min()[i]) * gen.uniform<double>() + box.min()[i];
  }
  return v;
}

void model::morph_into_saw_f77(rand_gen &gen,
                               const input::morph_into_saw_t &params) {

  auto min_dist_sq = pow(params.intersection_at, 2.0);

  using Vector = Eigen::Vector3d;

  double bond = 0.0;
  if (params.bond_distance.has_value()) {
    bond = params.bond_distance.value();
  } else {
    int num_native = 0;
    for (auto const &teth : tethers) {
      if (teth.length.has_value()) {
        num_native += 1;
        bond += teth.length.value();
      }
    }
    bond = bond / num_native;
  }

  Eigen::AlignedBox3d box;
  if (params.init_box_density > 0) {
    auto vol = (double)residues.size() / params.init_box_density;
    auto cell_a = cbrt(vol);
    box.min() = {-cell_a / 2.0, -cell_a / 2.0, -cell_a / 2.0};
    box.max() = -box.min();
    if (params.with_pbc) {
      Eigen::Vector3d box_dim = box.max() - box.min();
      model_box.set_cell(box_dim);
    }
  } else {
    if (params.with_pbc)
      throw std::runtime_error("with_pbc requires non-zero box");

    box.max() = box.min() = Eigen::Vector3d::Zero();
  }

  for (auto &chain : chains) {
    int n = (int)chain->residues.size();

    for (int kter = 0; kter < 9000; ++kter) {
      std::vector<double> phi(n), theta(n);

      auto pi = acos(-1.0);

      phi[0] = pi / 2.0;
      theta[0] = 0.0;

      phi[1] = 0.0;
      theta[1] = gen.uniform<double>() * pi / 3.0;

      for (int i = 2; i < n - 1; ++i) {
        phi[i] = (2.0 * gen.uniform<double>() - 1.0) * pi;
        theta[i] = gen.uniform<double>() * pi / 3.0;
      }

      std::vector<Eigen::Matrix3d> T(n);
      for (int i = 0; i < n - 1; ++i) {
        auto ct = cos(theta[i]), st = sin(theta[i]), cp = cos(phi[i]),
             sp = sin(phi[i]);

        T[i] << ct, st, 0.0, st * cp, -ct * cp, sp, st * sp, -ct * sp, -cp;
      }

      auto theta0 = acos(1.0 - 2.0 * gen.uniform<double>());
      auto phi0 = 2.0 * pi * gen.uniform<double>();

      std::vector<Vector> R(n);
      for (int i = 0; i < n - 1; ++i) {
        Vector r;
        r << bond * sin(theta0) * cos(phi0), bond * sin(theta0) * sin(phi0),
            bond * cos(theta0);

        for (int j = i; j >= 0; --j) {
          R[i] = T[j] * r;
          r = R[i];
        }
      }

      Vector ran = sample_from_box(box, gen);
      for (int i = 0; i < n; ++i) {
        Vector r = ran;
        for (int j = 0; j < i; ++j) {
          r += R[j];
        }
        chain->residues[i]->pos = r;
      }

      bool success = true;
      for (int i = 0; i < n - 3 && success; ++i) {
        auto r1 = chain->residues[i]->pos;
        for (int j = i + 3; j < n && success; ++j) {
          auto r2 = chain->residues[j]->pos;
          auto dx = !params.with_pbc ? (r2 - r1) : model_box.wrap(r1, r2);
          if (norm_squared(dx) < min_dist_sq) {
            success = false;
            break;
          }
        }
        if (!success)
          break;
      }

      if (success)
        break;
    }
  }

  if (!params.with_pbc) {
    box = Eigen::AlignedBox3d();
    for (auto const &res : residues)
      box.extend(convert<double>(res->pos));

    Vector ext = box.max() - box.min();
    model_box.set_cell(vec3d(ext));
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
