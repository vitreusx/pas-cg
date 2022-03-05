#include "langevin/legacy_step.h"
#include "utils/units.h"
using namespace cg::lang;

void legacy_step::operator()() const {
  auto local_gen = *gen;
  real noise_factor = (real)sqrt(2.0 * kB * temperature);
  solver_real dt_inv = 1.0 / dt;
  auto gamma_factor_sqrt = sqrt(gamma_factor);
  auto dt_sqrt = sqrt(dt);

  for (int idx = 0; idx < num_particles; ++idx) {
    auto aa_idx = (uint8_t)atype[idx];
    auto gamma = gamma_factor * mass[aa_idx];
    auto noise_sd = noise_factor * gamma_factor_sqrt * mass_rsqrt[aa_idx];

    auto [noise_x, noise_y] = local_gen.normal_x2<real>();
    auto noise_z = local_gen.normal<real>();
    auto noise = vec3r(noise_x, noise_y, noise_z);

    y1[idx] = y1[idx] + noise_sd * noise * dt * dt_sqrt;
    vec3sr a_ = F[idx] * mass_inv[aa_idx] - gamma * dt_inv * y1[idx];

    vec3sr error = y2[idx] - a_ * (dt * dt / 2.0);
    y0[idx] -= 3.0 / 16.0 * error;
    y1[idx] -= 251.0 / 360.0 * error;
    y2[idx] -= 1.0 * error;
    y3[idx] -= 11.0 / 18.0 * error;
    y4[idx] -= 1.0 / 6.0 * error;
    y5[idx] -= 1.0 / 60.0 * error;

    y0[idx] += y1[idx] + y2[idx] + y3[idx] + y4[idx] + y5[idx];
    y1[idx] += 2.0 * y2[idx] + 3.0 * y3[idx] + 4.0 * y4[idx] + 5.0 * y5[idx];
    y2[idx] += 3.0 * y3[idx] + 6.0 * y4[idx] + 10.0 * y5[idx];
    y3[idx] += 4.0 * y4[idx] + 10.0 * y5[idx];
    y4[idx] += 5.0 * y5[idx];

    r[idx] = y0[idx];
    v[idx] = y1[idx] * dt_inv;
  }

  *true_t += dt;
  *t = (real)*true_t;
  *gen = local_gen;
}