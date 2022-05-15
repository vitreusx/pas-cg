#pragma once

#ifndef COMPAT_MODE
#include "random/xorshift.h"
#else
#include "random/nr.h"
#endif
