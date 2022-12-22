#include <cg/langevin/fast_step.h>
#include <cg/utils/sanitize.h>
#include <cg/utils/units.h>

namespace cg::lang {
void fast_step::set_params(cg::real gamma_factor, cg::real temperature,
                           cg::solver_real dt) {
  this->dt = dt;

  dt_inv = (solver_real)1.0 / dt;

  gamma2 = gamma_factor * dt_inv;
  const2 = (real)2.0 * (real)kB * temperature * gamma_factor * dt;
  const2 = cg::sqrt(const2) * dt;
  deltsq = (dt * dt) / (solver_real)2.0;
}

void fast_step::omp_async() const {
  if (*t == 0) {
#pragma omp for schedule(static) nowait
    for (int idx = 0; idx < num_particles; ++idx)
      y2[idx] = F[idx] * deltsq;

#pragma omp barrier
  }

#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < num_particles; ++idx) {
    auto aa_idx = (uint8_t)atype[idx];
    auto gam2 = gamma2 * mass_inv[aa_idx];
    auto con2 = const2 * mass_rsqrt[aa_idx];

    vec3r f = F[idx] - gam2 * y1[idx];
    auto [noise_x, noise_y] = gens[idx].normal_x2<real>();
    auto noise_z = gens[idx].normal<real>();
    y1[idx] = y1[idx] + con2 * vec3r(noise_x, noise_y, noise_z);

    vec3sr err = y2[idx] - deltsq * f;
    y0[idx] -= f02 * err;
    y1[idx] -= f12 * err;
    y2[idx] -= 1.0 * err;
    y3[idx] -= f32 * err;
    y4[idx] -= f42 * err;
    y5[idx] -= f52 * err;

    y0[idx] += y1[idx] + y2[idx] + y3[idx] + y4[idx] + y5[idx];
    y1[idx] += 2.0 * y2[idx] + 3.0 * y3[idx] + 4.0 * y4[idx] + 5.0 * y5[idx];
    y2[idx] += 3.0 * y3[idx] + 6.0 * y4[idx] + 10.0 * y5[idx];
    y3[idx] += 4.0 * y4[idx] + 10.0 * y5[idx];
    y4[idx] += 5.0 * y5[idx];

    r[idx] = y0[idx];
    v[idx] = y1[idx] * dt_inv;
  }

#pragma omp single nowait
  {
    *true_t += dt;
    *t = (real)*true_t;
    ++*step_idx;
  }
}
} // namespace cg::lang