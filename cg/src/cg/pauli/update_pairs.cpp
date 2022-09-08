#include <cg/pauli/update_pairs.h>
namespace cg::pauli {

void update_pairs::operator()() const {
  pairs->clear();

  for (decltype(auto) nl_pair : nl->pairs) {
    if (nl_pair.taken())
      continue;

    auto i1 = nl_pair.i1(), i2 = nl_pair.i2();

    auto r1 = r[i1], r2 = r[i2];
    auto dist = norm(simul_box->wrap<vec3r>(r1, r2));
    if (dist < r_excl + nl->orig_pad) {
      pairs->emplace_back(i1, i2);
      nl_pair.taken() = true;
    }
  }
}
} // namespace cg::pauli