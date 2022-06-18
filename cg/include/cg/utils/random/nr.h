#pragma once
#include <algorithm>
#include <cg/types/vec3.h>
#include <cmath>
#include <cstdint>

namespace cg {
class rand_gen {
private:
  static constexpr int im1 = 2147483563, im2 = 2147483399, imm1 = im1 - 1,
                       ia1 = 40014, ia2 = 40692, iq1 = 53668, iq2 = 52774,
                       ir1 = 12211, ir2 = 3791, ntab = 32,
                       ndiv = 1 + imm1 / ntab;

  static constexpr double eps = 1.2e-7f, rnmx = 1.f - eps,
                          am = 1.f / (float)im1;

  int iy = 0, idum = -448, idum2 = 123456789;
  int iv[ntab];

public:
  rand_gen() : rand_gen(0) {}

  inline explicit rand_gen(int seed) {
    ran2(seed);
    ran2(seed);
  }

  inline double ran2(int iseed = 0) {
    int k, j;
    if (iseed != 0)
      idum = -iseed;

    if (idum <= 0) {
      idum2 = idum = std::max(-idum, 1);
      for (j = ntab + 7; j >= 0; --j) {
        k = idum / iq1;
        idum = ia1 * (idum - k * iq1) - k * ir1;
        if (idum < 0)
          idum += im1;
        if (j < ntab)
          iv[j] = idum;
      }
      iy = iv[0];
    }

    k = idum / iq1;
    idum = ia1 * (idum - k * iq1) - k * ir1;
    if (idum < 0)
      idum += im1;

    k = idum2 / iq2;
    idum2 = ia2 * (idum2 - k * iq2) - k * ir2;
    if (idum2 < 0)
      idum2 += im2;

    j = iy / ndiv;
    iy = iv[j] - idum2;
    iv[j] = idum;
    if (iy < 1)
      iy += imm1;

    return std::min(am * iy, rnmx);
  }

  inline uint64_t operator()() {
    double dval = ran2();
    return uint64_t(dval * (double)std::numeric_limits<uint64_t>::max());
  }

  inline rand_gen spawn() {
    uint64_t seed = rand_gen(*this)();
    seed = (seed ^ (seed >> 30)) * 0xBF58476D1CE4E5B9;
    seed = (seed ^ (seed >> 27)) * 0x94D049BB133111EB;
    return rand_gen(seed ^ (seed >> 31));
  }

  template <typename U> inline U uniform() {
    return ran2();
  }

  template <typename U> inline U uniform(U a, U b) {
    return (b - a) * uniform<U>() + a;
  }

  template <typename U> inline U normal() {
    auto r1 = uniform<U>(), r2 = uniform<U>();
    U r = sqrt((U)(-2.0) * log(r1));
    U t = (U)(2.0 * M_PI) * r2;
    return r * cos(t);
  }

  template <typename U> inline std::pair<U, U> normal_x2() {
    auto r1 = uniform<U>(), r2 = uniform<U>();
    U r = sqrt((U)(-2.0) * log(r1));
    U t = (U)(2.0 * M_PI) * r2;
    auto n1 = r * cos(t), n2 = r * sin(t);
    return std::make_pair(n1, n2);
  }

  template <typename U> inline vec3<U> sphere() {
    auto [x, y] = normal_x2<U>();
    auto z = normal<U>();
    return unit(vec3<U>(x, y, z));
  }
};
} // namespace cg