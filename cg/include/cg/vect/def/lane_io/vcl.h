#pragma once
#include "../lane/vcl.h"
#include "../repr.h"
#include "decl.h"
#include <immintrin.h>

namespace nitro::def {
// IMPL(int8_t, 8, 256, Vec8i);

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

// IMPL(uint8_t, 8, 256, Vec8ui);

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

// IMPL(int16_t, 8, 256, Vec8i);

inline void load(Vec8i &data, int16_t const *src) {
  Vec8s vec;
  vec.load(src);
  data = extend(vec);
}

inline void store(Vec8i const &data, int16_t *dst) {
  compress(data).store(dst);
}

// IMPL(uint16_t, 8, 256, Vec8ui);

inline void load(Vec8ui &data, uint16_t const *src) {
  Vec8us vec;
  vec.load(src);
  data = extend(vec);
}

inline void store(Vec8ui const &data, uint16_t *dst) {
  compress(data).store(dst);
}

// IMPL(int32_t, 8, 256, Vec8i);

inline void load(Vec8i &data, int32_t const *src) {
  data.load(src);
}

inline void store(Vec8i const &data, int32_t *dst) {
  data.store(dst);
}

// IMPL(uint32_t, 8, 256, Vec8ui);

inline void load(Vec8ui &data, uint32_t const *src) {
  data.load(src);
}

inline void store(Vec8ui const &data, uint32_t *dst) {
  data.store(dst);
}

// IMPL(float, 8, 256, Vec8f);

inline void load(Vec8f &data, float const *src) {
  data.load(src);
}

inline void store(Vec8f const &data, float *dst) {
  data.store(dst);
}

// IMPL(bool, 8, 256, Vec8ib);

inline void load(Vec8ib &data, bool const *src) {
  Vec8ui idata;
  load(idata, reinterpret_cast<repr_t<bool> const *>(src));
  data = (idata != 0);
}

inline void store(Vec8ib const &data, bool *dst) {
  static auto zero = Vec8ui(0), one = Vec8ui(1);
  Vec8ui idata = select(data, zero, one);
  store(idata, reinterpret_cast<repr_t<bool> *>(dst));
}

// IMPL(int8_t, 4, 256, Vec4q);

inline void load(Vec4q &data, int8_t const *src) {
  data = Vec4q(src[0], src[1], src[2], src[3]);
}

inline void store(Vec4q const &data, int8_t *dst) {
  int64_t values[4];
  data.store(values);
  for (int idx = 0; idx < 4; ++idx)
    dst[idx] = values[idx];
}

// IMPL(uint8_t, 4, 256, Vec4uq);

inline void load(Vec4uq &data, uint8_t const *src) {
  data = Vec4uq(src[0], src[1], src[2], src[3]);
}

inline void store(Vec4q const &data, uint8_t *dst) {
  uint64_t values[4];
  data.store(values);
  for (int idx = 0; idx < 4; ++idx)
    dst[idx] = values[idx];
}

// IMPL(int16_t, 4, 256, Vec4q);

inline void load(Vec4q &data, int16_t const *src) {
#ifdef __AVX2__
  __m128i small_load = _mm_loadl_epi64((const __m128i *)src);
  data = _mm256_cvtepi16_epi64(small_load);
#else
  data = Vec4q(src[0], src[1], src[2], src[3]);
#endif
}

inline void store(Vec4q const &data, int16_t *dst) {
  int64_t values[4];
  data.store(values);
  for (int idx = 0; idx < 4; ++idx)
    dst[idx] = values[idx];
}

// IMPL(uint16_t, 4, 256, Vec4uq);

inline void load(Vec4uq &data, uint16_t const *src) {
#ifdef __AVX2__
  __m128i small_load = _mm_loadl_epi64((const __m128i *)src);
  data = _mm256_cvtepu16_epi64(small_load);
#else
  data = Vec4q(src[0], src[1], src[2], src[3]);
#endif
}

inline void store(Vec4uq const &data, uint16_t *dst) {
  uint64_t values[4];
  data.store(values);
  for (int idx = 0; idx < 4; ++idx)
    dst[idx] = values[idx];
}

// IMPL(int32_t, 4, 256, Vec4q);

inline void load(Vec4q &data, int32_t const *src) {
  Vec4i vec;
  vec.load(src);
  data = extend(vec);
}

inline void store(Vec4q const &data, int32_t *dst) {
  compress(data).store(dst);
}

// IMPL(uint32_t, 4, 256, Vec4uq);

inline void load(Vec4uq &data, uint32_t const *src) {
  Vec4ui vec;
  vec.load(src);
  data = extend(vec);
}

inline void store(Vec4uq const &data, uint32_t *dst) {
  compress(data).store(dst);
}

// IMPL(int64_t, 4, 256, Vec4q);

inline void load(Vec4q &data, int64_t const *src) {
  data.load(src);
}

inline void store(Vec4q const &data, int64_t *dst) {
  data.store(dst);
}

// IMPL(uint64_t, 4, 256, Vec4uq);

inline void load(Vec4uq &data, uint64_t const *src) {
  data.load(src);
}

inline void store(Vec4uq const &data, uint64_t *dst) {
  data.store(dst);
}

// IMPL(float, 4, 256, Vec4d);

inline void load(Vec4d &data, float const *src) {
  Vec4f vec;
  vec.load(src);
  data = to_double(vec);
}

inline void store(Vec4d const &data, float *dst) {
  to_float(data).store(dst);
}

// IMPL(double, 4, 256, Vec4d);

inline void load(Vec4d &data, double const *src) {
  data.load(src);
}

inline void store(Vec4d const &data, double *dst) {
  data.store(dst);
}

// IMPL(bool, 4, 256, Vec4ib);

inline void load(Vec4qb &data, bool const *src) {
  Vec4uq idata;
  load(idata, reinterpret_cast<repr_t<bool> const *>(src));
  data = (idata != 0);
}

inline void store(Vec4qb const &data, bool *dst) {
  Vec4uq idata = select(data, Vec4uq(0), Vec4uq(1));
  store(idata, reinterpret_cast<repr_t<bool> *>(dst));
}

// IMPL(int8_t, 16, 512, Vec16i);

inline void load(Vec16i &data, int8_t const *src) {
  Vec16c vec;
  vec.load(src);
  data = extend(extend(vec));
}

inline void store(Vec16i const &data, int8_t *dst) {
  compress(compress(data)).store(dst);
}

// IMPL(uint8_t, 16, 512, Vec16ui);

inline void load(Vec16ui &data, uint8_t const *src) {
  Vec16uc vec;
  vec.load(src);
  data = extend(extend(vec));
}

inline void store(Vec16ui const &data, uint8_t *dst) {
  compress(compress(data)).store(dst);
}

// IMPL(int16_t, 16, 512, Vec16i);

inline void load(Vec16i &data, int16_t const *src) {
  Vec16s vec;
  vec.load(src);
  data = extend(vec);
}

inline void store(Vec16i const &data, int16_t *dst) {
  compress(data).store(dst);
}

// IMPL(uint16_t, 16, 512, Vec16ui);

inline void load(Vec16ui &data, uint16_t const *src) {
  Vec16us vec;
  vec.load(src);
  data = extend(vec);
}

inline void store(Vec16ui const &data, uint16_t *dst) {
  compress(data).store(dst);
}

// IMPL(int32_t, 16, 512, Vec16i);

inline void load(Vec16i &data, int32_t const *src) {
  data.load(src);
}

inline void store(Vec16i const &data, int32_t *dst) {
  data.store(dst);
}

// IMPL(uint32_t, 16, 512, Vec16ui);

inline void load(Vec16ui &data, uint32_t const *src) {
  data.load(src);
}

inline void store(Vec16ui const &data, uint32_t *dst) {
  data.store(dst);
}

// IMPL(float, 16, 512, Vec16f);

inline void load(Vec16f &data, float const *src) {
  data.load(src);
}

inline void store(Vec16f const &data, float *dst) {
  data.store(dst);
}

// IMPL(bool, 16, 512, Vec16ib);

inline void load(Vec16ib &data, bool const *src) {
  Vec16ui idata;
  load(idata, reinterpret_cast<repr_t<bool> const *>(src));
  data = (idata != 0);
}

inline void store(Vec16ib const &data, bool *dst) {
  Vec16ui idata = select(data, Vec16ui(0), Vec16ui(1));
  store(idata, reinterpret_cast<repr_t<bool> *>(dst));
}

// IMPL(int8_t, 8, 512, Vec8q);

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

// IMPL(uint8_t, 8, 512, Vec8uq);

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

// IMPL(int16_t, 8, 512, Vec8q);

inline void load(Vec8q &data, int16_t const *src) {
  Vec8s vec;
  vec.load(src);
  data = extend(extend(vec));
}

inline void store(Vec8q const &data, int16_t *dst) {
  compress(compress(data)).store(dst);
}

// IMPL(uint16_t, 8, 512, Vec8uq);

inline void load(Vec8uq &data, uint16_t const *src) {
  Vec8us vec;
  vec.load(src);
  data = extend(extend(vec));
}

inline void store(Vec8uq &data, uint16_t *dst) {
  compress(compress(data)).store(dst);
}

// IMPL(int32_t, 8, 512, Vec8q);

inline void load(Vec8q &data, int32_t const *src) {
  Vec8i vec;
  vec.load(src);
  data = extend(vec);
}

inline void store(Vec8q const &data, int32_t *dst) {
  compress(data).store(dst);
}

// IMPL(uint32_t, 8, 512, Vec8uq);

inline void load(Vec8uq &data, uint32_t const *src) {
  Vec8ui vec;
  vec.load(src);
  data = extend(vec);
}

inline void store(Vec8uq &data, uint32_t *dst) {
  compress(data).store(dst);
}

// IMPL(int64_t, 8, 512, Vec8q);

inline void load(Vec8q &data, int64_t const *src) {
  data.load(src);
}

inline void store(Vec8q const &data, int64_t *dst) {
  data.store(dst);
}

// IMPL(uint64_t, 8, 512, Vec8uq);

inline void load(Vec8uq &data, uint64_t const *src) {
  data.load(src);
}

inline void store(Vec8uq const &data, uint64_t *dst) {
  data.store(dst);
}

// IMPL(float, 8, 512, Vec8d);

inline void load(Vec8d &data, float const *src) {
  Vec8f proxy;
  proxy.load(src);
  data = to_double(proxy);
}

inline void store(Vec8d const &data, float *dst) {
  to_float(data).store(dst);
}

// IMPL(double, 8, 512, Vec8d);

inline void load(Vec8d &data, double const *src) {
  data.load(src);
}

inline void store(Vec8d const &data, double *dst) {
  data.store(dst);
}

// IMPL(bool, 8, 512, Vec8qb);

inline void load(Vec8qb &data, bool const *src) {
  Vec8uq idata;
  load(idata, reinterpret_cast<repr_t<bool> const *>(src));
  data = (idata != 0);
}

inline void store(Vec8qb const &data, bool *dst) {
  Vec8uq idata = select(data, Vec8uq(0), Vec8uq(1));
  store(idata, reinterpret_cast<repr_t<bool> *>(dst));
}

template <typename Lane, typename T,
          typename = std::enable_if_t<is_vcl_lane_v<Lane>>>
Lane _construct(selector<Lane>, T const *src) {
  Lane lane;
  load(lane, src);
  return lane;
}

template <typename Lane, typename T, typename Idx, std::size_t... I>
Lane _gather_aux(selector<Lane>, T const *src, Idx const *idxes,
                 std::index_sequence<I...>) {
  return Lane(src[idxes[I]]...);
}

template <typename Lane, typename T, typename Idxes,
          typename = std::enable_if_t<is_vcl_lane_v<Lane>>>
Lane _gather(selector<Lane>, T const *src, Idxes const &idxes) {
  lane_type_t<Idxes> idxes_[lane_size_v<Idxes>];
  store(idxes, idxes_);
  return _gather_aux(selector<Lane>{}, src, idxes_,
                     std::make_index_sequence<lane_size_v<Lane>>{});
}

template <typename Lane, typename T, typename Idxes,
          typename = std::enable_if_t<is_vcl_lane_v<Lane>>>
void scatter(Lane const &data, T *dst, Idxes const &idxes) {
  lane_type_t<Idxes> idxes_[lane_size_v<Idxes>];
  store(idxes, idxes_);
  for (int idx = 0; idx < lane_size_v<Lane>; ++idx)
    dst[idxes_[idx]] = data[idx];
}

} // namespace nitro::def