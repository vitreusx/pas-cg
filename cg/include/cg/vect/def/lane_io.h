#pragma once
#include "lane.h"

namespace nitro::def {
inline void load(Vec8ui &data, uint8_t const *src) {
#ifdef __AVX2__
  __m128i small_load = _mm_loadl_epi64((const __m128i *)src);
  data = _mm256_cvtepu8_epi32(small_load);
#else
  data = Vec8ui(src[0], src[1], src[2], src[3], src[4], src[5], src[6], src[7]);
#endif
}

inline void store(Vec8ui const &data, uint8_t *dst) {
  _mm_storel_epi64((__m128i *)dst, compress(data));
}

inline void load(Vec8ui &data, uint16_t const *src) {
  Vec8us vec;
  vec.load(src);
  data = extend(vec);
}

inline void store(Vec8ui const &data, uint16_t *dst) {
  compress(data).store(dst);
}

inline void load(Vec8ui &data, uint32_t const *src) {
  data.load(src);
}

inline void store(Vec8ui const &data, uint32_t *dst) {
  data.store(dst);
}

inline void load(Vec8uq &data, uint8_t const *src) {
#ifdef __AVX512F__
  __m128i small_load = _mm_loadl_epi64((const __m128i *)ptr);
  data = _mm512_cvtepu8_epi64(small_load);
#else
  Vec8ui proxy;
  proxy.load(src);
  data = extend(proxy);
#endif
}

inline void store(Vec8uq const &data, uint8_t *dst) {
  _mm_storel_epi64((__m128i *)dst, compress(compress(data)));
}

inline void load(Vec8uq &data, uint16_t const *src) {
  Vec8us vec;
  vec.load(src);
  data = extend(extend(vec));
}

inline void store(Vec8uq &data, uint16_t *dst) {
  compress(compress(data)).store(dst);
}

inline void load(Vec8uq &data, uint32_t const *src) {
  Vec8ui vec;
  vec.load(src);
  data = extend(vec);
}

inline void store(Vec8uq &data, uint32_t *dst) {
  compress(data).store(dst);
}

inline void load(Vec8uq &data, uint64_t const *src) {
  data.load(src);
}

inline void store(Vec8uq const &data, uint64_t *dst) {
  data.store(dst);
}

inline void load(Vec8i &data, int8_t const *src) {
#ifdef __AVX2__
  __m128i small_load = _mm_loadl_epi64((const __m128i *)src);
  data = _mm256_cvtepi8_epi32(small_load);
#else
  data = Vec8i(src[0], src[1], src[2], src[3], src[4], src[5], src[6], src[7]);
#endif
}

inline void store(Vec8i const &data, int8_t *dst) {
  _mm_storel_epi64((__m128i *)dst, compress(data));
}

inline void load(Vec8i &data, int16_t const *src) {
  Vec8s vec;
  vec.load(src);
  data = extend(vec);
}

inline void store(Vec8i const &data, int16_t *dst) {
  compress(data).store(dst);
}

inline void load(Vec8i &data, int32_t const *src) {
  data.load(src);
}

inline void store(Vec8i const &data, int32_t *dst) {
  data.store(dst);
}

inline void load(Vec8q &data, int8_t const *src) {
#ifdef __AVX512F__
  __m128i small_load = _mm_loadl_epi64((const __m128i *)src);
  data = _mm512_cvtepi8_epi64(small_load);
#else
  Vec8i proxy;
  load(proxy, src);
  data = extend(proxy);
#endif
}

inline void store(Vec8q &data, int8_t *dst) {
  _mm_storel_epi64((__m128i *)dst, compress(compress(data)));
}

inline void load(Vec8q &data, int16_t const *src) {
  Vec8s vec;
  vec.load(src);
  data = extend(extend(vec));
}

inline void store(Vec8q const &data, int16_t *dst) {
  compress(compress(data)).store(dst);
}

inline void load(Vec8q &data, int32_t const *src) {
  Vec8i vec;
  vec.load(src);
  data = extend(vec);
}

inline void store(Vec8q const &data, int32_t *dst) {
  compress(data).store(dst);
}

inline void load(Vec8q &data, int64_t const *src) {
  data.load(src);
}

inline void store(Vec8q const &data, int64_t *dst) {
  data.store(dst);
}

inline void load(Vec16ui &data, uint8_t const *src) {
  Vec16uc vec;
  vec.load(src);
  data = extend(extend(vec));
}

inline void store(Vec16ui const &data, uint8_t *dst) {
  compress(compress(data)).store(dst);
}

inline void load(Vec16ui &data, uint16_t const *src) {
  Vec16us vec;
  vec.load(src);
  data = extend(vec);
}

inline void store(Vec16ui const &data, uint16_t *dst) {
  compress(data).store(dst);
}

inline void load(Vec16ui &data, uint32_t const *src) {
  data.load(src);
}

inline void store(Vec16ui const &data, uint32_t *dst) {
  data.store(dst);
}

inline void load(Vec16i &data, int8_t const *src) {
  Vec16c vec;
  vec.load(src);
  data = extend(extend(vec));
}

inline void store(Vec16i const &data, int8_t *dst) {
  compress(compress(data)).store(dst);
}

inline void load(Vec16i &data, int16_t const *src) {
  Vec16s vec;
  vec.load(src);
  data = extend(vec);
}

inline void store(Vec16i const &data, int16_t *dst) {
  compress(data).store(dst);
}

inline void load(Vec16i &data, int32_t const *src) {
  data.load(src);
}

inline void store(Vec16i const &data, int32_t *dst) {
  data.store(dst);
}

inline void load(Vec8f &data, float const *src) {
  data.load(src);
}

inline void store(Vec8f const &data, float *dst) {
  data.store(dst);
}

inline void load(Vec8d &data, float const *src) {
  Vec8f proxy;
  proxy.load(src);
  data = to_double(proxy);
}

inline void store(Vec8d const &data, float *dst) {
  to_float(data).store(dst);
}

inline void load(Vec8d &data, double const *src) {
  data.load(src);
}

inline void store(Vec8d const &data, double *dst) {
  data.store(dst);
}

inline void load(Vec16f &data, float const *src) {
  data.load(src);
}

inline void store(Vec16f const &data, float *dst) {
  data.store(dst);
}

inline void load(Vec8ib &data, bool const *src) {
  Vec8ui idata;
  load(idata, reinterpret_cast<uint_repr<bool> const *>(src));
  data = (idata != 0);
}

inline void store(Vec8ib const &data, bool *dst) {
  static auto zero = Vec8ui(0), one = Vec8ui(1);
  Vec8ui idata = select(data, zero, one);
  store(idata, reinterpret_cast<uint_repr<bool> *>(dst));
}

inline void load(Vec8qb &data, bool const *src) {
  Vec8uq idata;
  load(idata, reinterpret_cast<uint_repr<bool> const *>(src));
  data = (idata != 0);
}

inline void store(Vec8qb const &data, bool *dst) {
  Vec8uq idata = select(data, Vec8q(0), Vec8q(1));
  store(idata, reinterpret_cast<uint_repr<bool> *>(dst));
}

inline void load(Vec16ib &data, bool const *src) {
  Vec16ui idata;
  load(idata, reinterpret_cast<uint_repr<bool> const *>(src));
  data = (idata != 0);
}

inline void store(Vec16ib const &data, bool *dst) {
  Vec16ui idata = select(data, Vec16ui(0), Vec16ui(1));
  store(idata, reinterpret_cast<uint_repr<bool> *>(dst));
}

template <typename L, typename T>
void load(L &data, T const *src) {
  load(data, reinterpret_cast<uint_repr<T> const *>(src));
}

template <typename L, typename T>
void store(L const &data, T *dst) {
  store(data, reinterpret_cast<uint_repr<T> *>(dst));
}
} // namespace nitro::def