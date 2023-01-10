#include <cg/dh/const_cuda.h>

namespace cg::const_dh {
void eval_forces_cuda::set_V_factor(real factor) {
  V_factor = 1.0 / (4.0 * M_PI * factor);
}

__global__ static void kernel(eval_forces_cuda eval) {
  auto idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= eval.es_pairs.size())
    return;

  auto es = eval.es_pairs[idx];
  auto i1 = es.i1(), i2 = es.i2();
  auto q1_x_q2 = es.q1_x_q2();

  auto r1 = eval.r[i1], r2 = eval.r[i2];
  auto r12 = eval.simul_box.wrap<vec3r>(r1, r2);

  auto r12_n = norm(r12);
  if (r12_n > eval.cutoff)
    return;
  auto r12_rn = (real)1.0f / r12_n;
  auto r12_u = r12 * r12_rn;

  auto Vij =
      eval.V_factor * q1_x_q2 * exp(-r12_n * eval.screen_dist_inv) * r12_rn;
  auto dVij_dr = -Vij * (eval.screen_dist_inv + r12_rn);

  atomicAdd(eval.V, Vij);

  auto f = r12_u * dVij_dr;
  eval.F[i1].cudaAtomicAdd(f);
  eval.F[i2].cudaAtomicSub(f);
}

void eval_forces_cuda::operator()(int block_size, cudaStream_t stream) const {
  auto n = es_pairs.size();
  dim3 block(block_size), grid(n / block_size + 1);
  kernel<<<grid, block, 0, stream>>>(*this);
}

} // namespace cg::const_dh