#include "features/langevin/step.h"
#include "utils/units.h"
using namespace cg::lang;

void step::operator()() const {
  auto local_gen = *gen;
  real noise_factor = (real)sqrt(2.0 * kB * temperature);
  solver_real dt_inv = 1.0 / dt;
  auto gamma_factor_sqrt = sqrt(gamma_factor);

  for (int idx = 0; idx < num_particles; ++idx) {
    auto gamma = gamma_factor * mass->at(idx);
    auto noise_sd = noise_factor * gamma_factor_sqrt * mass_rsqrt->at(idx);

    auto [noise_x, noise_y] = local_gen.normal_x2<real>();
    auto noise_z = local_gen.normal<real>();
    auto noise = vec3r(noise_x, noise_y, noise_z);

    auto a_ =
        F->at(idx) * mass_inv->at(idx) - gamma * v->at(idx) + noise_sd * noise;

    vec3sr error = y2->at(idx) - a_ * (dt * dt / 2.0);
    y0->at(idx) -= 3.0 / 16.0 * error;
    y1->at(idx) -= 251.0 / 360.0 * error;
    y2->at(idx) -= 1.0 * error;
    y3->at(idx) -= 11.0 / 18.0 * error;
    y4->at(idx) -= 1.0 / 6.0 * error;
    y5->at(idx) -= 1.0 / 60.0 * error;

    y0->at(idx) +=
        y1->at(idx) + y2->at(idx) + y3->at(idx) + y4->at(idx) + y5->at(idx);
    y1->at(idx) += 2.0 * y2->at(idx) + 3.0 * y3->at(idx) + 4.0 * y4->at(idx) +
                   5.0 * y5->at(idx);
    y2->at(idx) += 3.0 * y3->at(idx) + 6.0 * y4->at(idx) + 10.0 * y5->at(idx);
    y3->at(idx) += 4.0 * y4->at(idx) + 10.0 * y5->at(idx);
    y4->at(idx) += 5.0 * y5->at(idx);

    r->at(idx) = y0->at(idx);
    v->at(idx) = y1->at(idx) * dt_inv;
  }

  *true_t += dt;
  *t = (real)*true_t;
  *gen = local_gen;
}