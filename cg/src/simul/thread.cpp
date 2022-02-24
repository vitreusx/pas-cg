#include "simul/thread_state.h"
using namespace cg::simul;

thread_state::thread_state(int num_residues, uint64_t seed)
    : dyn(num_residues), gen(seed){};