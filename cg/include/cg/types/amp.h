#pragma once
#include <cg/types/vec3.h>
#include <cstddef>

namespace cg {
#ifndef COMPAT_MODE
using real = float;
#else
using real = double;
#endif

using solver_real = double;

using vec3r = vec3<real>;
using vec3sr = vec3<solver_real>;
} // namespace cg