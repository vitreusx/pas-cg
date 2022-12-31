#include <cg/nl/cell_update.h>
#include <cg/utils/sort.h>
#include <limits>

namespace cg::nl {
template <typename E1, typename E2>
static bool operator<(nl::pair_expr<E1> const &left,
                      nl::pair_expr<E2> const &right) {
  return std::make_pair(left.i1(), left.i2()) <
         std::make_pair(right.i1(), right.i2());
}

void cell_update::omp_async() const {
  auto inf = std::numeric_limits<real>::infinity();
  vec3r r_min(inf, inf, inf), r_max(-inf, -inf, -inf);

#pragma omp single nowait
  {
    *shared_r_min = vec3r(inf, inf, inf);
    *shared_r_max = vec3r(-inf, -inf, -inf);
  };

#pragma omp for schedule(static)
  for (auto idx : idxes) {
    auto r_ = simul_box->wrap((vec3r)r[idx]);

    r_min.x() = min(r_min.x(), r_.x());
    r_min.y() = min(r_min.y(), r_.y());
    r_min.z() = min(r_min.z(), r_.z());

    r_max.x() = max(r_max.x(), r_.x());
    r_max.y() = max(r_max.y(), r_.y());
    r_max.z() = max(r_max.z(), r_.z());
  }

#pragma omp critical
  {
    shared_r_min->x() = min(shared_r_min->x(), r_min.x());
    shared_r_min->y() = min(shared_r_min->y(), r_min.y());
    shared_r_min->z() = min(shared_r_min->z(), r_min.z());

    shared_r_max->x() = max(shared_r_max->x(), r_max.x());
    shared_r_max->y() = max(shared_r_max->y(), r_max.y());
    shared_r_max->z() = max(shared_r_max->z(), r_max.z());
  };

#pragma omp barrier

  r_min = *shared_r_min;
  r_max = *shared_r_max;

  vec3r ext = r_max - r_min;

  auto eff_cutoff = cutoff + pad;
  auto cell_a = eff_cutoff;

  auto nx = max((int)std::floor(ext.x() / cell_a), 1);
  auto ny = max((int)std::floor(ext.y() / cell_a), 1);
  auto nz = max((int)std::floor(ext.z() / cell_a), 1);

  auto cell_x = ext.x() / nx, cell_x_inv = (real)1.0 / cell_x;
  auto cell_y = ext.y() / ny, cell_y_inv = (real)1.0 / cell_y;
  auto cell_z = ext.z() / nz, cell_z_inv = (real)1.0 / cell_z;

  auto max_dx = (int)std::ceil(eff_cutoff * cell_x_inv);
  max_dx = clamp(max_dx, 0, nx - 1);
  auto max_dy = (int)std::ceil(eff_cutoff * cell_y_inv);
  max_dy = clamp(max_dy, 0, ny - 1);
  auto max_dz = (int)std::ceil(eff_cutoff * cell_z_inv);
  max_dz = clamp(max_dz, 0, nz - 1);

  offsets.clear();
  auto eff_cutoff_sq = eff_cutoff * eff_cutoff;
  for (int dz = -max_dz; dz <= max_dz; ++dz) {
    auto diff_z = cell_z * max((int)abs(dz) - 1, 0);
    for (int dy = -max_dy; dy <= max_dy; ++dy) {
      auto diff_y = cell_y * max((int)abs(dy) - 1, 0);
      for (int dx = -max_dx; dx <= max_dx; ++dx) {
        auto diff_x = cell_x * max((int)abs(dx) - 1, 0);
        auto diff_v = vec3r(diff_x, diff_y, diff_z);

        if (norm_squared(diff_v) <= eff_cutoff_sq)
          offsets.emplace_back(dx, dy, dz);
      }
    }
  }

#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < idxes.size(); ++idx) {
    auto r_ = simul_box->wrap((vec3r)r[idxes[idx]]);
    auto ix = (int)std::floor((r_.x() - r_min.x()) * cell_x_inv);
    ix = clamp(ix, 0, nx - 1);
    auto iy = (int)std::floor((r_.y() - r_min.y()) * cell_y_inv);
    iy = clamp(iy, 0, ny - 1);
    auto iz = (int)std::floor((r_.z() - r_min.z()) * cell_z_inv);
    iz = clamp(iz, 0, nz - 1);
    auto cell_idx = ix + nx * (iy + ny * iz);
    cell_idxes[idx] = std::make_pair(cell_idx, idx);
  }

#pragma omp barrier

#pragma omp single nowait
  {
    std::sort(cell_idxes.begin(), cell_idxes.end());

    unique_cells->clear();
    unique_cell_idxes->clear();

    auto cur_idx = -1;
    int *end_ptr = nullptr;
    for (int idx = 0; idx < cell_idxes.size(); ++idx) {
      auto cell_idx = cell_idxes[idx].first;
      if (cell_idx != cur_idx) {
        if (end_ptr)
          *end_ptr = idx;

        auto &[cell_idx_, beg, end] =
            unique_cells->emplace_back(cell_idx, idx, -1);
        end_ptr = &end;
        cur_idx = cell_idx;

        (*unique_cell_idxes)[cell_idx] = unique_cells->size() - 1;
      }
    }

    if (end_ptr)
      *end_ptr = cell_idxes.size();

    *total_pairs = *cur_pairs_offset = 0;
    nl_data->pairs.clear();
  }

#pragma omp barrier

  local_pairs.clear();

#pragma omp for schedule(static) nowait
  for (auto const &[idx1, beg1, end1] : *unique_cells) {
    auto ix1 = idx1 % nx, iy1 = (idx1 / nx) % ny, iz1 = (idx1 / nx) / ny;
    for (auto [dx, dy, dz] : offsets) {
      auto ix2 = ix1 + dx, iy2 = iy1 + dy, iz2 = iz1 + dz;

      if (simul_box->cell_inv.x() == (real)0.0) {
        if (ix2 < 0 || ix2 >= nx)
          continue;
      } else {
        ix2 = (ix2 + nx) % nx;
      }

      if (simul_box->cell_inv.y() == (real)0.0) {
        if (iy2 < 0 || iy2 >= ny)
          continue;
      } else {
        iy2 = (iy2 + ny) % ny;
      }

      if (simul_box->cell_inv.z() == (real)0.0) {
        if (iz2 < 0 || iz2 >= nz)
          continue;
      } else {
        iz2 = (iz2 + nz) % nz;
      }

      auto idx2 = ix2 + nx * (iy2 + ny * iz2);
      if (idx1 > idx2)
        continue;

      auto iter = unique_cell_idxes->find(idx2);
      if (iter == unique_cell_idxes->end())
        continue;

      auto [_, beg2, end2] = unique_cells->at(iter->second);
      for (int i1 = beg1; i1 < end1; ++i1) {
        auto res_idx1 = idxes[cell_idxes[i1].second];
        auto r1 = r[res_idx1];
        auto chain1 = chain_idx[res_idx1], seq1 = seq_idx[res_idx1];

        auto beg2_ = (idx1 != idx2) ? beg2 : (i1 + 1);
        for (int i2 = beg2_; i2 < end2; ++i2) {
          auto res_idx2 = idxes[cell_idxes[i2].second];
          auto r2 = r[res_idx2];
          auto chain2 = chain_idx[res_idx2], seq2 = seq_idx[res_idx2];

          auto diff = seq1 > seq2 ? seq1 - seq2 : seq2 - seq1;
          if (chain1 == chain2 && diff < 3)
            continue;

          auto orig_dist = norm(simul_box->wrap((vec3r)(r2 - r1)));
          if (orig_dist < eff_cutoff) {
            auto res1_ = res_idx1, res2_ = res_idx2;
            if (res1_ > res2_)
              std::swap(res1_, res2_);
            local_pairs.emplace_back(res1_, res2_, orig_dist, false);
          }
        }
      }
    }
  }

#pragma omp critical
  *total_pairs += local_pairs.size();

#pragma omp barrier

#pragma omp single nowait
  {
    auto dummy = nl::pair(0, 0, (real)0.0, false);
    nl_data->pairs.resize(*total_pairs, dummy);
  };

#pragma omp barrier

  int offset = 0;
#pragma omp critical
  {
    offset = *cur_pairs_offset;
    *cur_pairs_offset += local_pairs.size();
  };

  for (int idx = 0; idx < local_pairs.size(); ++idx) {
    nl_data->pairs[idx + offset] = local_pairs[idx];
  }

#pragma omp barrier

  sort_omp_async(nl_data->pairs.begin(), nl_data->pairs.end());
#pragma omp barrier

#pragma omp single nowait
  {
    //    std::sort(nl_data->pairs.begin(), nl_data->pairs.end());

    nl_data->orig_pad = pad;
    for (int idx = 0; idx < r.size(); ++idx)
      nl_data->orig_r[idx] = r[idx];
    nl_data->orig_pbc = *simul_box;
    *invalid = false;
  };
}
} // namespace cg::nl