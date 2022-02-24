#pragma once
#include <cg/amino/amino_acid.h>

namespace cg::qa {
class contact_type {
public:
  inline contact_type() : val{decltype(val)(NONE())} {};

  static inline contact_type NONE() { return contact_type(0); }

  static inline contact_type BACK_BACK() { return contact_type(1); }

  static inline contact_type BACK_SIDE() { return contact_type(2); }

  static inline contact_type SIDE_BACK() { return contact_type(3); }

  static inline contact_type SIDE_SIDE(amino_acid const &a1,
                                       amino_acid const &a2) {

    int16_t val = (int16_t)4u + (int16_t)(uint8_t)a1 +
                  (int16_t)(uint8_t)a2 * (int16_t)amino_acid::NUM_TYPES;
    return contact_type(val);
  }

  explicit inline operator int16_t() const { return val; }

  static const int NUM_TYPES =
      4 + amino_acid::NUM_TYPES * amino_acid::NUM_TYPES;

private:
  explicit inline contact_type(int16_t val) : val{val} {};

  int16_t val;
};
} // namespace cg::qa