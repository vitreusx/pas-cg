#include "qa/contact_type.h"
using namespace cg::qa;

bool contact_type::operator==(const contact_type &other) const {
  return val == other.val;
}

bool contact_type::operator!=(const contact_type &other) const {
  return val != other.val;
}
