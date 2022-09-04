#include <cg/base_forces/shifted_lj.h>
#include <cg/qa/loop_over_candidates.h>

namespace cg::qa {

void loop_over_candidates::operator()() const {
  for (int idx = 0; idx < free_pairs->size(); ++idx) {
    iter(idx);
  }
}

void loop_over_candidates::iter(int idx) const {
  decltype(auto) node = free_pairs->at(idx);

  if (node.is_vacant())
    return;

  auto pair = node.item();
  auto orig_dist = pair.orig_dist();
  if (orig_dist - *total_disp > max_req_dist)
    return;

  auto i1 = pair.i1(), i2 = pair.i2();

  auto r1 = r[i1], r2 = r[i2];
  auto r12 = simul_box->wrap(r1, r2);
  auto r12_rn = norm_inv(r12);
  auto r12_u = r12 * r12_rn;

  if ((real)1.0 > r12_rn * max_req_dist)
    return;

  if (1.0 < r12_rn * rep_cutoff) {
    auto [V_, dV_dr] = shifted_lj(rep_depth, rep_cutoff)(r12_rn);
    *V += V_;
    F[i1] += r12_u * dV_dr;
    F[i2] -= r12_u * dV_dr;
  }

  auto n1 = n[i1], n2 = n[i2], h1 = h[i1], h2 = h[i2];

  auto back1 = (abs(dot(h1, r12_u)) > min_abs_cos_hr);
  auto back2 = (abs(dot(h2, r12_u)) > min_abs_cos_hr);
  auto bb_extra = (abs(dot(h1, h2)) > min_abs_cos_hh);
  auto side1 = (dot(n1, r12_u) < max_cos_nr);
  auto side2 = (-dot(n2, r12_u) < max_cos_nr);

  auto atype1 = atype[i1], atype2 = atype[i2];

  contact_type type;
  if (back1 && back2 && bb_extra)
    type = contact_type::BACK_BACK();
  else if (back1 && side2)
    type = contact_type::BACK_SIDE();
  else if (side1 && back2)
    type = contact_type::SIDE_BACK();
  else if (side1 && side2)
    type = contact_type::SIDE_SIDE(atype1, atype2);

  if ((short)type == (short)contact_type::NONE())
    return;

  auto req_dist = req_min_dist[(short)type] * ((real)1.0 + formation_tolerance);
  if (1.0f > r12_rn * req_dist)
    return;

  auto ptype2 = ptype[(uint8_t)atype2];
  sync_data sync_diff1;
  sync_diff1.back() += back1 ? 1 : 0;
  sync_diff1.side_all() += side1 ? 1 : 0;
  sync_diff1.side_polar() += (side1 && (ptype2 == POLAR)) ? 1 : 0;
  sync_diff1.side_hydrophobic() += (side1 && (ptype2 == HYDROPHOBIC)) ? 1 : 0;

  sync_data sync1_after_formation = sync[i1] - sync_diff1;
  if (!sync1_after_formation.is_valid())
    return;

  auto ptype1 = ptype[(uint8_t)atype1];
  sync_data sync_diff2;
  sync_diff2.back() += back2 ? 1 : 0;
  sync_diff2.side_all() += side2 ? 1 : 0;
  sync_diff2.side_polar() += (side2 && (ptype1 == POLAR)) ? 1 : 0;
  sync_diff2.side_hydrophobic() += (side2 && (ptype1 == HYDROPHOBIC)) ? 1 : 0;

  sync_data sync2_after_formation = sync[i2] - sync_diff2;
  if (!sync2_after_formation.is_valid())
    return;

  static auto ss_type = contact_type::SIDE_SIDE(aa_code::CYS, aa_code::CYS);
  if (disulfide_special_criteria && (int16_t)type == (int16_t)ss_type) {
    if (neigh[i1] + neigh[i2] > max_neigh_count || part_of_ssbond[i1] ||
        part_of_ssbond[i2])
      return;
  }

#pragma omp critical
  candidates->emplace_back(i1, i2, orig_dist, (real)1.0 / r12_rn, idx, type,
                           sync_diff1, sync_diff2);
}

void loop_over_candidates::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < free_pairs->size(); ++idx) {
    iter(idx);
  }
}

void loop_over_candidates::for_slice(int from, int to) const {
  for (int idx = from; idx < to; ++idx)
    iter(idx);
}

int loop_over_candidates::total_size() const {
  return free_pairs->size();
}

} // namespace cg::qa