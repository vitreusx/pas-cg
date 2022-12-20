#pragma once
#include "../lane/vcl.h"
#include "../repr.h"
#include "decl.h"
#include <immintrin.h>
#include <iostream>

namespace nitro::def {
// IMPL(int8_t, 8, 256, Vec8i);

inline void vcl_load(Vec8i &data, int8_t const *src) {
#ifdef __AVX2__
  __m128i small_load = _mm_loadl_epi64((const __m128i *)src);
  data = _mm256_cvtepi8_epi32(small_load);
#else
  data = Vec8i(src[0], src[1], src[2], src[3], src[4], src[5], src[6], src[7]);
#endif
}

inline void vcl_store(Vec8i const &data, int8_t *dst) {
  _mm_storel_epi64((__m128i *)dst, compress(data));
}

// IMPL(uint8_t, 8, 256, Vec8ui);

inline void vcl_load(Vec8ui &data, uint8_t const *src) {
#ifdef __AVX2__
  __m128i small_load = _mm_loadl_epi64((const __m128i *)src);
  data = _mm256_cvtepu8_epi32(small_load);
#else
  data = Vec8ui(src[0], src[1], src[2], src[3], src[4], src[5], src[6], src[7]);
#endif
}

inline void vcl_store(Vec8ui const &data, uint8_t *dst) {
  _mm_storel_epi64((__m128i *)dst, compress(data));
}

// IMPL(int16_t, 8, 256, Vec8i);

inline void vcl_load(Vec8i &data, int16_t const *src) {
  Vec8s vec;
  vec.load(src);
  data = extend(vec);
}

inline void vcl_store(Vec8i const &data, int16_t *dst) {
  compress(data).store(dst);
}

// IMPL(uint16_t, 8, 256, Vec8ui);

inline void vcl_load(Vec8ui &data, uint16_t const *src) {
  Vec8us vec;
  vec.load(src);
  data = extend(vec);
}

inline void vcl_store(Vec8ui const &data, uint16_t *dst) {
  compress(data).store(dst);
}

// IMPL(int32_t, 8, 256, Vec8i);

inline void vcl_load(Vec8i &data, int32_t const *src) {
  data.load(src);
}

inline void vcl_store(Vec8i const &data, int32_t *dst) {
  data.store(dst);
}

// IMPL(uint32_t, 8, 256, Vec8ui);

inline void vcl_load(Vec8ui &data, uint32_t const *src) {
  data.load(src);
}

inline void vcl_store(Vec8ui const &data, uint32_t *dst) {
  data.store(dst);
}

// IMPL(float, 8, 256, Vec8f);

inline void vcl_load(Vec8f &data, float const *src) {
  data.load(src);
}

inline void vcl_store(Vec8f const &data, float *dst) {
  data.store(dst);
}

// IMPL(bool, 8, 256, Vec8ib);

inline void vcl_load(Vec8ib &data, bool const *src) {
  // Vec8ui idata;
  // load(idata, reinterpret_cast<repr_t<bool> const *>(src));
  // data = (idata != 0);
  data = Vec8ib(src[0], src[1], src[2], src[3], src[4], src[5], src[6], src[7]);
}

inline void vcl_store(Vec8ib const &data, bool *dst) {
  // Vec8ui idata = select(data, Vec8ui(1), Vec8ui(0));
  // store(idata, reinterpret_cast<repr_t<bool> *>(dst));
  uint32_t idata[8];
  vcl_store(select(data, Vec8ui(1), Vec8ui(0)), idata);
  for (int idx = 0; idx < 8; ++idx)
    dst[idx] = idata[idx];
}

// IMPL(int8_t, 4, 256, Vec4q);

inline void vcl_load(Vec4q &data, int8_t const *src) {
  data = Vec4q(src[0], src[1], src[2], src[3]);
}

inline void vcl_store(Vec4q const &data, int8_t *dst) {
  int64_t values[4];
  data.store(values);
  for (int idx = 0; idx < 4; ++idx)
    dst[idx] = values[idx];
}

// IMPL(uint8_t, 4, 256, Vec4uq);

inline void vcl_load(Vec4uq &data, uint8_t const *src) {
  data = Vec4uq(src[0], src[1], src[2], src[3]);
}

inline void vcl_store(Vec4q const &data, uint8_t *dst) {
  uint64_t values[4];
  data.store(values);
  for (int idx = 0; idx < 4; ++idx)
    dst[idx] = values[idx];
}

// IMPL(int16_t, 4, 256, Vec4q);

inline void vcl_load(Vec4q &data, int16_t const *src) {
#ifdef __AVX2__
  __m128i small_load = _mm_loadl_epi64((const __m128i *)src);
  data = _mm256_cvtepi16_epi64(small_load);
#else
  data = Vec4q(src[0], src[1], src[2], src[3]);
#endif
}

inline void vcl_store(Vec4q const &data, int16_t *dst) {
  int64_t values[4];
  data.store(values);
  for (int idx = 0; idx < 4; ++idx)
    dst[idx] = values[idx];
}

// IMPL(uint16_t, 4, 256, Vec4uq);

inline void vcl_load(Vec4uq &data, uint16_t const *src) {
#ifdef __AVX2__
  __m128i small_load = _mm_loadl_epi64((const __m128i *)src);
  data = _mm256_cvtepu16_epi64(small_load);
#else
  data = Vec4q(src[0], src[1], src[2], src[3]);
#endif
}

inline void vcl_store(Vec4uq const &data, uint16_t *dst) {
  uint64_t values[4];
  data.store(values);
  for (int idx = 0; idx < 4; ++idx)
    dst[idx] = values[idx];
}

// IMPL(int32_t, 4, 256, Vec4q);

inline void vcl_load(Vec4q &data, int32_t const *src) {
  Vec4i vec;
  vec.load(src);
  data = extend(vec);
}

inline void vcl_store(Vec4q const &data, int32_t *dst) {
  compress(data).store(dst);
}

// IMPL(uint32_t, 4, 256, Vec4uq);

inline void vcl_load(Vec4uq &data, uint32_t const *src) {
  Vec4ui vec;
  vec.load(src);
  data = extend(vec);
}

inline void vcl_store(Vec4uq const &data, uint32_t *dst) {
  compress(data).store(dst);
}

// IMPL(int64_t, 4, 256, Vec4q);

inline void vcl_load(Vec4q &data, int64_t const *src) {
  data.load(src);
}

inline void vcl_store(Vec4q const &data, int64_t *dst) {
  data.store(dst);
}

// IMPL(uint64_t, 4, 256, Vec4uq);

inline void vcl_load(Vec4uq &data, uint64_t const *src) {
  data.load(src);
}

inline void vcl_store(Vec4uq const &data, uint64_t *dst) {
  data.store(dst);
}

// IMPL(float, 4, 256, Vec4d);

inline void vcl_load(Vec4d &data, float const *src) {
  Vec4f vec;
  vec.load(src);
  data = to_double(vec);
}

inline void vcl_store(Vec4d const &data, float *dst) {
  to_float(data).store(dst);
}

// IMPL(double, 4, 256, Vec4d);

inline void vcl_load(Vec4d &data, double const *src) {
  data.load(src);
}

inline void store(Vec4d const &data, double *dst) {
  data.store(dst);
}

// IMPL(bool, 4, 256, Vec4ib);

inline void vcl_load(Vec4qb &data, bool const *src) {
  // Vec4uq idata;
  // load(idata, reinterpret_cast<repr_t<bool> const *>(src));
  // data = (idata != 0);
  data = Vec4qb(src[0], src[1], src[2], src[3]);
}

inline void vcl_store(Vec4qb const &data, bool *dst) {
  // Vec4uq idata = select(data, Vec4uq(0), Vec4uq(1));
  // store(idata, reinterpret_cast<repr_t<bool> *>(dst));
  uint64_t idata[4];
  vcl_store(select(data, Vec4uq(1), Vec4uq(0)), idata);
  for (int idx = 0; idx < 4; ++idx)
    dst[idx] = idata[idx];
}

// IMPL(int8_t, 16, 512, Vec16i);

inline void vcl_load(Vec16i &data, int8_t const *src) {
  Vec16c vec;
  vec.load(src);
  data = extend(extend(vec));
}

inline void vcl_store(Vec16i const &data, int8_t *dst) {
  compress(compress(data)).store(dst);
}

// IMPL(uint8_t, 16, 512, Vec16ui);

inline void vcl_load(Vec16ui &data, uint8_t const *src) {
  Vec16uc vec;
  vec.load(src);
  data = extend(extend(vec));
}

inline void vcl_store(Vec16ui const &data, uint8_t *dst) {
  compress(compress(data)).store(dst);
}

// IMPL(int16_t, 16, 512, Vec16i);

inline void vcl_load(Vec16i &data, int16_t const *src) {
  Vec16s vec;
  vec.load(src);
  data = extend(vec);
}

inline void vcl_store(Vec16i const &data, int16_t *dst) {
  compress(data).store(dst);
}

// IMPL(uint16_t, 16, 512, Vec16ui);

inline void vcl_load(Vec16ui &data, uint16_t const *src) {
  Vec16us vec;
  vec.load(src);
  data = extend(vec);
}

inline void vcl_store(Vec16ui const &data, uint16_t *dst) {
  compress(data).store(dst);
}

// IMPL(int32_t, 16, 512, Vec16i);

inline void vcl_load(Vec16i &data, int32_t const *src) {
  data.load(src);
}

inline void vcl_store(Vec16i const &data, int32_t *dst) {
  data.store(dst);
}

// IMPL(uint32_t, 16, 512, Vec16ui);

inline void vcl_load(Vec16ui &data, uint32_t const *src) {
  data.load(src);
}

inline void vcl_store(Vec16ui const &data, uint32_t *dst) {
  data.store(dst);
}

// IMPL(float, 16, 512, Vec16f);

inline void vcl_load(Vec16f &data, float const *src) {
  data.load(src);
}

inline void vcl_store(Vec16f const &data, float *dst) {
  data.store(dst);
}

// IMPL(bool, 16, 512, Vec16ib);

inline void vcl_load(Vec16ib &data, bool const *src) {
  // Vec16ui idata;
  // load(idata, reinterpret_cast<repr_t<bool> const *>(src));
  // data = (idata != 0);
  data = Vec16ib(src[0], src[1], src[2], src[3], src[4], src[5], src[6], src[7],
                 src[8], src[9], src[10], src[11], src[12], src[13], src[14],
                 src[15]);
}

inline void vcl_store(Vec16ib const &data, bool *dst) {
  // Vec16ui idata = select(data, Vec16ui(0), Vec16ui(1));
  // store(idata, reinterpret_cast<repr_t<bool> *>(dst));
  uint32_t idata[16];
  vcl_store(select(data, Vec16ui(1), Vec16ui(0)), idata);
  for (int idx = 0; idx < 16; ++idx)
    dst[idx] = idata[idx];
}

// IMPL(int8_t, 8, 512, Vec8q);

inline void vcl_load(Vec8q &data, int8_t const *src) {
#ifdef __AVX512F__
  __m128i small_load = _mm_loadl_epi64((const __m128i *)src);
  data = _mm512_cvtepi8_epi64(small_load);
#else
  Vec8i proxy;
  vcl_load(proxy, src);
  data = extend(proxy);
#endif
}

inline void vcl_store(Vec8q &data, int8_t *dst) {
  _mm_storel_epi64((__m128i *)dst, compress(compress(data)));
}

// IMPL(uint8_t, 8, 512, Vec8uq);

inline void vcl_load(Vec8uq &data, uint8_t const *src) {
#ifdef __AVX512F__
  __m128i small_load = _mm_loadl_epi64((const __m128i *)src);
  data = _mm512_cvtepu8_epi64(small_load);
#else
  Vec8ui proxy;
  proxy.load(src);
  data = extend(proxy);
#endif
}

inline void vcl_store(Vec8uq const &data, uint8_t *dst) {
  _mm_storel_epi64((__m128i *)dst, compress(compress(data)));
}

// IMPL(int16_t, 8, 512, Vec8q);

inline void vcl_load(Vec8q &data, int16_t const *src) {
  Vec8s vec;
  vec.load(src);
  data = extend(extend(vec));
}

inline void vcl_store(Vec8q const &data, int16_t *dst) {
  compress(compress(data)).store(dst);
}

// IMPL(uint16_t, 8, 512, Vec8uq);

inline void vcl_load(Vec8uq &data, uint16_t const *src) {
  Vec8us vec;
  vec.load(src);
  data = extend(extend(vec));
}

inline void vcl_store(Vec8uq &data, uint16_t *dst) {
  compress(compress(data)).store(dst);
}

// IMPL(int32_t, 8, 512, Vec8q);

inline void vcl_load(Vec8q &data, int32_t const *src) {
  Vec8i vec;
  vec.load(src);
  data = extend(vec);
}

inline void vcl_store(Vec8q const &data, int32_t *dst) {
  compress(data).store(dst);
}

// IMPL(uint32_t, 8, 512, Vec8uq);

inline void vcl_load(Vec8uq &data, uint32_t const *src) {
  Vec8ui vec;
  vec.load(src);
  data = extend(vec);
}

inline void vcl_store(Vec8uq &data, uint32_t *dst) {
  compress(data).store(dst);
}

// IMPL(int64_t, 8, 512, Vec8q);

inline void vcl_load(Vec8q &data, int64_t const *src) {
  data.load(src);
}

inline void vcl_store(Vec8q const &data, int64_t *dst) {
  data.store(dst);
}

// IMPL(uint64_t, 8, 512, Vec8uq);

inline void vcl_load(Vec8uq &data, uint64_t const *src) {
  data.load(src);
}

inline void vcl_store(Vec8uq const &data, uint64_t *dst) {
  data.store(dst);
}

// IMPL(float, 8, 512, Vec8d);

inline void vcl_load(Vec8d &data, float const *src) {
  Vec8f proxy;
  proxy.load(src);
  data = to_double(proxy);
}

inline void vcl_store(Vec8d const &data, float *dst) {
  to_float(data).store(dst);
}

// IMPL(double, 8, 512, Vec8d);

inline void vcl_load(Vec8d &data, double const *src) {
  data.load(src);
}

inline void vcl_store(Vec8d const &data, double *dst) {
  data.store(dst);
}

// IMPL(bool, 8, 512, Vec8qb);
#if INSTRSET < 10
inline void vcl_load(Vec8qb &data, bool const *src) {
  // Vec8uq idata;
  // load(idata, reinterpret_cast<repr_t<bool> const *>(src));
  // data = (idata != 0);
  data = Vec8qb(src[0], src[1], src[2], src[3], src[4], src[5], src[6], src[7]);
}

inline void vcl_store(Vec8qb const &data, bool *dst) {
  // Vec8uq idata = select(data, Vec8uq(0), Vec8uq(1));
  // store(idata, reinterpret_cast<repr_t<bool> *>(dst));
  uint64_t idata[8];
  vcl_store(select(data, Vec8uq(1), Vec8uq(0)), idata);
  for (int idx = 0; idx < 8; ++idx)
    dst[idx] = idata[idx];
}
#endif

template <typename Lane, typename T>
Lane vcl_construct(T const *src) {
  Lane lane;
  vcl_load(lane, src);
  return lane;
}

template <typename Lane, typename T, typename Idx, std::size_t... I>
Lane vcl_gather_aux(T const *src, Idx const *idxes, std::index_sequence<I...>) {
  return Lane(src[idxes[I]]...);
}

template <typename Lane, typename T, typename Idxes,
          typename = std::enable_if_t<is_vcl_lane_v<Lane>>>
Lane vcl_gather(T const *src, Idxes const &idxes) {
  lane_type_t<Idxes> idxes_[lane_size_v<Idxes>];
  vcl_store(idxes, idxes_);
  return vcl_gather_aux<Lane>(src, idxes_,
                              std::make_index_sequence<lane_size_v<Lane>>{});
  //  return Lane(lookup<lane_size_v<Lane>>(idxes, src));
}

template <typename Lane, typename T, typename Idx, typename Mask,
          std::size_t... I>
Lane vcl_masked_gather_aux(T const *src, Idx const *idxes, Mask const *mask,
                           std::index_sequence<I...>) {
  return Lane((mask[I] ? src[idxes[I]] : T())...);
}

template <typename Lane, typename T, typename Idxes, typename Mask,
          typename = std::enable_if_t<is_vcl_lane_v<Lane>>>
Lane vcl_masked_gather(T const *src, Idxes const &idxes, Mask const &mask) {
  lane_type_t<Idxes> idxes_[lane_size_v<Idxes>];
  vcl_store(idxes, idxes_);
  lane_type_t<Mask> mask_[lane_size_v<Mask>];
  vcl_store(mask, mask_);

  return vcl_masked_gather_aux<Lane>(
      src, idxes_, mask_, std::make_index_sequence<lane_size_v<Lane>>{});
}

template <typename Lane, typename T, typename Idxes,
          typename = std::enable_if_t<is_vcl_lane_v<Lane>>>
void vcl_scatter(Lane const &data, T *dst, Idxes const &idxes) {
  lane_type_t<Idxes> idxes_[lane_size_v<Idxes>];
  vcl_store(idxes, idxes_);
  for (size_t idx = 0; idx < lane_size_v<Lane>; ++idx)
    dst[idxes_[idx]] = data[idx];
  //  scatter(idxes, (uint32_t)0xffffffff, data, dst);
}

template <typename Lane, typename T, typename Idxes, typename Mask,
          typename = std::enable_if_t<is_vcl_lane_v<Lane>>>
void vcl_masked_scatter(Lane const &data, T *dst, Idxes const &idxes,
                        Mask const &mask) {
  lane_type_t<Idxes> idxes_[lane_size_v<Idxes>];
  vcl_store(idxes, idxes_);
  lane_type_t<Mask> mask_[lane_size_v<Mask>];
  vcl_store(mask, mask_);
  for (size_t idx = 0; idx < lane_size_v<Idxes>; ++idx)
    if (mask_[idx])
      dst[idxes_[idx]] = data[idx];
}

} // namespace nitro::def