#pragma once

#undef MIXED_PRECISION
#define MIXED_PRECISION true

#undef LEGACY_MODE
#define LEGACY_MODE true

#ifndef MIXED_PRECISION
#define MIXED_PRECISION false
#endif

#ifndef LEGACY_MODE
#define LEGACY_MODE false
#endif