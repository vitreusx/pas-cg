#pragma once
#include "data.h"
#include "dynamics.h"
#include <cg/base_forces/sink_lj.h>
#include <cg/cuda/dev_var.h>
#include <cg/cuda/stream.h>

namespace cg::simul {
class cuda_state {
public:
  vect::cuda_vector<vec3r> r, F;
  cuda::dev_var<real> V;
  vect::cuda_vector<dh::pair> es_pairs;
  vect::cuda_vector<pid::bundle> bundles;
  vect::cuda_vector<int> prev, next;
  vect::cuda_vector<amino_acid> atype;
  vect::cuda_vector<sink_lj> pid_ss_ljs;

public:
  dynamics staging;
  cuda::stream stream;

public:
  void reset_dyn(int block_size, cudaStream_t stream = 0);
  void pull_dyn(int block_size, cudaStream_t stream = 0);
};
} // namespace cg::simul