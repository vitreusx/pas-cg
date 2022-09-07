#include <Eigen/Geometry>
#include <cg/input/model.h>
#include <cg/sbox/box.h>
#include <cg/sbox/pbc.h>
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

  cryst1 = other.cryst1;

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

static Eigen::Vector3d sample_from_box(Eigen::AlignedBox3d const &box,
                                       rand_gen &gen) {
  Eigen::Vector3d v;
  for (int i = 0; i < 3; ++i) {
    v[i] = (box.max()[i] - box.min()[i]) * gen.uniform<double>() + box.min()[i];
  }
  return v;
}

void model::morph_into_saw(rand_gen &gen,
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

  sbox::pbc<double> saw_pbc;

  Eigen::AlignedBox3d start_box;
  if (std::holds_alternative<input::morph_into_saw_t::start_box_params>(
          params.start_box)) {
    auto const &box_params =
        std::get<input::morph_into_saw_t::start_box_params>(params.start_box);

    double cell_a = 0.0;
    if (box_params.size.has_value()) {
      cell_a = box_params.size.value();
    } else if (box_params.density.has_value()) {
      auto density = box_params.density.value();
      auto vol = (double)residues.size() / density;
      cell_a = cbrt(vol);
    }

    start_box.min() = {-cell_a / 2.0, -cell_a / 2.0, -cell_a / 2.0};
    start_box.max() = -start_box.min();

    if (params.pbc) {
      Eigen::Vector3d box_dim = start_box.max() - start_box.min();
      saw_pbc.set_cell(box_dim);

      saw_box = sbox::box<double>();
      saw_box->min = start_box.min();
      saw_box->max = start_box.max();
    }
  } else if (std::holds_alternative<input::morph_into_saw_t::start_box_origin>(
                 params.start_box)) {
    start_box.min() = start_box.max() = Eigen::Vector3d::Zero();
    if (params.pbc)
      throw std::runtime_error("SAW conformation generation with periodic "
                               "boundary conditions requires non-zero box");
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

      Vector ran = sample_from_box(start_box, gen);
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
          auto dx = saw_pbc.wrap<vec3<double>>(r1, r2);
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

void model::remove_native_structure() {
  for (auto &ang : angles)
    ang.theta.reset();
  for (auto &dih : dihedrals)
    dih.phi.reset();
  contacts = {};
}
} // namespace cg::input
