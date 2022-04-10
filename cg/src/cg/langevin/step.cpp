#include <cg/langevin/step.h>
#include <cg/utils/math.h>
#include <cg/utils/sanitize.h>
#include <cg/utils/units.h>
namespace cg::lang {

void step::operator()() const {
  auto local_gen = *gen;
  real noise_factor = (real)sqrt(2.0 * kB * temperature);
  solver_real dt_inv = 1.0 / dt;
  auto gamma_factor_sqrt = sqrt(gamma_factor);

  for (int idx = 0; idx < num_particles; ++idx) {
    auto aa_idx = (uint8_t)atype[idx];
    auto gamma = gamma_factor * mass[aa_idx];
    auto noise_sd = noise_factor * gamma_factor_sqrt * mass_rsqrt[aa_idx];

    auto [noise_x, noise_y] = local_gen.normal_x2<real>();
    auto noise_z = local_gen.normal<real>();
    auto noise = vec3r(noise_x, noise_y, noise_z);

    vec3r f = F[idx];
    sanitize(f, (real)1e3);
    auto a_ = f * mass_inv[aa_idx] - gamma * v[idx] + noise_sd * noise;

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

void step::omp_async() const {
  auto local_gen = *gen;
  real noise_factor = (real)sqrt(2.0 * kB * temperature);
  solver_real dt_inv = 1.0 / dt;
  auto gamma_factor_sqrt = sqrt(gamma_factor);

#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < num_particles; ++idx) {
    auto aa_idx = (uint8_t)atype[idx];
    auto gamma = gamma_factor * mass[aa_idx];
    auto noise_sd = noise_factor * gamma_factor_sqrt * mass_rsqrt[aa_idx];

    auto [noise_x, noise_y] = local_gen.normal_x2<real>();
    auto noise_z = local_gen.normal<real>();
    auto noise = vec3r(noise_x, noise_y, noise_z);

    vec3r f = F[idx];
    sanitize(f, (real)1e3);
    auto a_ = f * mass_inv[aa_idx] - gamma * v[idx] + noise_sd * noise;

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

#pragma omp master
  {
    *true_t += dt;
    *t = (real)*true_t;
  }
  *gen = local_gen;
}
} // namespace cg::lang