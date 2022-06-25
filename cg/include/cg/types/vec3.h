#pragma once

#include "vec3_def.h"

#ifdef NDEBUG
#include "vec3_ops_expr.h"
#else
#include "vec3_ops_debug.h"
#endif