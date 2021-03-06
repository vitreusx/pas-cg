#include <algorithm>
#include <cg/nl/cell_update.h>
#include <cg/nl/order.h>
#include <cg/utils/math.h>
namespace cg::nl {

using vec3i = cg::vec3<int>;

static int flat_idx(vec3i idx, vec3i max) {
  return idx.x() + max.x() * (idx.y() + max.y() * idx.z());
}

void cell_update::operator()() const {
  real min_real = std::numeric_limits<real>::lowest();
  real max_real = std::numeric_limits<real>::max();
  vec3r min_r(max_real, max_real, max_real);
  vec3r max_r(min_real, min_real, min_real);

  for (int idx = 0; idx < r.size(); ++idx) {
    auto pos = simul_box->wrap(r[idx]);
    min_r.x() = min(min_r.x(), pos.x());
    max_r.x() = max(max_r.x(), pos.x());
    min_r.y() = min(min_r.y(), pos.y());
    max_r.y() = max(max_r.y(), pos.y());
    min_r.z() = min(min_r.z(), pos.z());
    max_r.z() = max(max_r.z(), pos.z());
  }

  auto extent = max_r - min_r;
  auto radius = cutoff + pad;

  vec3r apx_num_cells = extent / radius;
  vec3i num_cells = {(int)floor(apx_num_cells.x()),
                     (int)floor(apx_num_cells.y()),
                     (int)floor(apx_num_cells.z())};

  vec3r cell_dim = {extent.x() / (real)num_cells.x(),
                    extent.y() / (real)num_cells.y(),
                    extent.z() / (real)num_cells.z()};

  vec3r cell_dim_inv = {(real)1.0 / cell_dim.x(), (real)1.0 / cell_dim.y(),
                        (real)1.0 / cell_dim.z()};

  auto total_num_cells = num_cells.x() * num_cells.y() * num_cells.z();
  num_res_in_cell->resize(total_num_cells);
  cell_offset->resize(total_num_cells);

  for (int idx = 0; idx < total_num_cells; ++idx) {
    num_res_in_cell->at(idx) = 0;
  }

  for (int idx = 0; idx < r.size(); ++idx) {
    auto rel_pos = simul_box->wrap(r[idx]) - min_r;
    vec3i res_coords = {(int)floor(rel_pos.x() * cell_dim_inv.x()),
                        (int)floor(rel_pos.y() * cell_dim_inv.y()),
                        (int)floor(rel_pos.z() * cell_dim_inv.z())};

    res_coords = {clamp(res_coords.x(), 0, num_cells.x() - 1),
                  clamp(res_coords.y(), 0, num_cells.y() - 1),
                  clamp(res_coords.z(), 0, num_cells.z() - 1)};

    auto cur_res_cell_idx = flat_idx(res_coords, num_cells);
    res_cell_idx[idx] = cur_res_cell_idx;
    ++num_res_in_cell->at(cur_res_cell_idx);
  }

  int cur_cell_offset = 0;
  for (int idx = 0; idx < total_num_cells; ++idx) {
    cell_offset->at(idx) = cur_cell_offset;
    cur_cell_offset += num_res_in_cell->at(idx);
    num_res_in_cell->at(idx) = 0;
  }

  for (int idx = 0; idx < r.size(); ++idx) {
    auto cell_idx = res_cell_idx[idx];
    auto cur_reord_idx =
        cell_offset->at(cell_idx) + num_res_in_cell->at(cell_idx);
    reordered_idx[cur_reord_idx] = idx;
    ++num_res_in_cell->at(cell_idx);
  }

  all_pairs->clear();
  real radius_inv = (real)1.0 / radius;
  for (int x0 = 0; x0 < num_cells.x(); ++x0) {
    for (int y0 = 0; y0 < num_cells.y(); ++y0) {
      for (int z0 = 0; z0 < num_cells.z(); ++z0) {
        auto flat0 = flat_idx(vec3i(x0, y0, z0), num_cells);
        for (int dx = -1; dx <= 1; ++dx) {
          auto x1 = (x0 + dx + num_cells.x()) % num_cells.x();
          for (int dy = -1; dy <= 1; ++dy) {
            auto y1 = (y0 + dy + num_cells.y()) % num_cells.y();
            for (int dz = -1; dz <= 1; ++dz) {
              auto z1 = (z0 + dz + num_cells.z()) % num_cells.z();
              auto flat1 = flat_idx(vec3i(x1, y1, z1), num_cells);

              for (int re0 = cell_begin(flat0); re0 < cell_end(flat0); ++re0) {
                auto i0 = reordered_idx[re0];
                auto r0 = r[i0];
                auto chain0 = chain_idx[i0], seq0 = seq_idx[i0];

                for (int re1 = cell_begin(flat1); re1 < cell_end(flat1);
                     ++re1) {
                  auto i1 = reordered_idx[re1];
                  if (i1 >= i0)
                    continue;

                  auto r1 = r[i1];
                  auto chain1 = chain_idx[i1], seq1 = seq_idx[i1];

                  auto diff = seq0 > seq1 ? seq0 - seq1 : seq1 - seq0;
                  if (chain0 == chain1 && diff < 3)
                    continue;

                  auto orig_dist_inv = norm_inv(simul_box->wrap(r0, r1));
                  if (orig_dist_inv < radius_inv)
                    continue;

                  all_pairs->emplace_back(i0, i1, (real)1.0 / orig_dist_inv);
                }
              }
            }
          }
        }
      }
    }
  }

  std::sort(
      all_pairs->begin(), all_pairs->end(), [](auto p, auto q) -> auto{
        return std::make_pair(p.i1(), p.i2()) < std::make_pair(q.i1(), q.i2());
      });

  nl_data->native.clear();
  nl_data->non_native.clear();

  int nat_cont_idx = 0;

  for (int idx = 0; idx < all_pairs->size(); ++idx) {
    bool non_native = true;
    auto cur_pair = all_pairs->at(idx);

    while (nat_cont_idx < all_nat_cont.size()) {
      auto cur_nat_cont = all_nat_cont[nat_cont_idx];

      if (cur_nat_cont < cur_pair) {
        ++nat_cont_idx;
      } else {
        if (!(cur_nat_cont > cur_pair)) {
          nl_data->native.push_back(cur_pair);
          non_native = false;
        }
        break;
      }
    }

    if (non_native) {
      nl_data->non_native.push_back(cur_pair);
    }
  }

  nl_data->orig_pad = pad;
  for (int idx = 0; idx < r.size(); ++idx)
    nl_data->orig_r[idx] = r[idx];
  nl_data->orig_pbc = *simul_box;
  *invalid = false;
}

int cell_update::cell_begin(int cell_idx) const {
  return cell_offset->at(cell_idx);
}

int cell_update::cell_end(int cell_idx) const {
  return cell_offset->at(cell_idx) + num_res_in_cell->at(cell_idx);
}
} // namespace cg::nl