#include <cg/qa/contact_type.h>
namespace cg::qa {

bool contact_type::operator==(const contact_type &other) const {
  return val == other.val;
}

bool contact_type::operator!=(const contact_type &other) const {
  return val != other.val;
}

std::vector<contact_type> contact_type::all() {
  std::vector<contact_type> ctypes(NUM_TYPES);
  for (int16_t ct_idx = 0; ct_idx < NUM_TYPES; ++ct_idx)
    ctypes[ct_idx] = contact_type(ct_idx);
  return ctypes;
}
} // namespace cg::qa