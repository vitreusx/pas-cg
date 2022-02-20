
/**
 * Force the Startup classes to run
 * at static initialization time:
 */
#define STATIC_CPARSE_STARTUP
#include "./builtin-features.inc"

#include <iostream>

struct foo {
    foo() {
        std::cout << "What the fuck" << '\n';
    }
} __foo;