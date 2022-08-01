#include <cg/wall/log_forces.h>

namespace cg::wall {
void log_forces::operator()() const {
  for (auto *wall : walls) {
    if (!wall)
      continue;

    wall->avg_F.add(*t, *wall->F);
    wall->work += dot(*wall->F, wall->shift);
  }

  auto neg_z_force = dot(*neg_z->F, neg_z->plane.normal());
  auto pos_z_force = dot(*pos_z->F, pos_z->plane.normal());
  auto z_force = neg_z_force + pos_z_force;
  avg_z_force->add(*t, z_force);
}
} // namespace cg::wall