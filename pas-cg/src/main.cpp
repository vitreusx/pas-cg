#include <iostream>
#include <pas-cg/utils/vec3.h>
using namespace cg;

template <typename E>
std::ostream &operator<<(std::ostream &os, vec3_expr<E> const &e) {
  os << "(" << e.x() << ", " << e.y() << ", " << e.z() << ")";
  return os;
}

int main() {
  std::cout << cross(vec3f(0, 1, 2), vec3f(3, 4, 5)) << '\n';
  return 0;
}