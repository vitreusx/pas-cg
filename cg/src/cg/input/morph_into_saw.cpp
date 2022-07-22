#include <cg/input/morph_into_saw.h>

namespace cg::input {

void morph_into_saw_t::link(ioxx::xyaml::proxy &p) {
  p["perform"] >> perform;

  auto box_p = p["start box"];
  if (box_p.IsScalar() && box_p.Scalar() == "origin") {
    start_box.emplace<start_box_origin>();
  } else {
    start_box_params params;
    box_p["density"] >> params.density;
    box_p["size"] >> params.size;
    start_box.emplace<start_box_params>(params);
  }

  p["intersection at"] >> intersection_at;
  p["num of retries"] >> num_of_retries;
  p["bond distance"] >> bond_distance;
  p["periodic boundary conditions"] >> pbc;
}
} // namespace cg::input