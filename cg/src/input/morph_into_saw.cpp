#include "input/morph_into_saw.h"
#include "utils/ioxx_interop.h"
namespace cg::input {

void morph_into_saw_t::link(ioxx::xyaml::proxy &p) {
  p["perform"] & perform;
  p["bond distance"] & bond_distance;
  p["residue density"] & residue_density;
  p["intersection at"] & intersection_at;
  p["num of retries"] & num_of_retries;
  p["infer simulation box"] & infer_box;
}
} // namespace cg::input