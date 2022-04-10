#include <cg/qa/contact_type.h>
namespace cg::qa {

bool contact_type::operator==(const contact_type &other) const {
  return val == other.val;
}

bool contact_type::operator!=(const contact_type &other) const {
  return val != other.val;
}
} // namespace cg::qa