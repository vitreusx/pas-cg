#pragma once
#include "decl.h"
#include <memory>

namespace nitro {
template <typename T> using def_allocator = std::allocator<T>;
}