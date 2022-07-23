#include <cg/wall/log_forces.h>

namespace cg::wall {
void log_forces::operator()() const {
  auto neg_z_perp_force = dot(*neg_z_force, neg_z_plane->normal());
  auto pos_z_perp_force = dot(*pos_z_force, pos_z_plane->normal());
  auto z_force = neg_z_perp_force + pos_z_perp_force;
  avg_z_force->add(*t, z_force);
}
} // namespace cg::wall