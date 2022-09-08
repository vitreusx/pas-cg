#include <cg/dh/update_pairs.h>
namespace cg::dh {

void update_pairs::operator()() const {
  pairs->clear();

  auto min_norm_inv = (real)1.0 / (cutoff + nl->orig_pad);

  for (decltype(auto) nl_pair : nl->pairs) {
    if (nl_pair.taken())
      continue;

    auto i1 = nl_pair.i1(), i2 = nl_pair.i2();
    auto q1_x_q2 = q[(uint8_t)atype[i1]] * q[(uint8_t)atype[i2]];
    if (q1_x_q2 == 0)
      continue;

    auto r1 = r[i1], r2 = r[i2];

    if (norm_inv(simul_box->wrap<vec3r>(r1, r2)) > min_norm_inv) {
      pairs->emplace_back(i1, i2, q1_x_q2);
      //      nl_pair.taken() = true;
    }
  }
}
} // namespace cg::dh