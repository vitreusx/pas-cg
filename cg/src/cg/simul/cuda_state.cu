#include <cg/simul/cuda_state.h>

namespace cg::simul {
static __global__ void reset_F_kernel(vect::view<vec3r> F, int n) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= n)
    return;

  F[idx] = vec3r::Zero();
}

void cuda_state::reset_dyn(int block_size, cudaStream_t stream) {
  V.set((real)0, stream);

  int n = F.size();
  dim3 grid(n / block_size + 1), block(block_size);
  reset_F_kernel<<<grid, block, 0, stream>>>(F, n);
}

void cuda_state::pull_dyn(int block_size, cudaStream_t stream) {
  V.fetch(stream);
  staging.V = (real)V;
  F.push_to(staging.F, stream);
}
} // namespace cg::simul