#pragma once
#include <algorithm>
#include <iterator>
#include <omp.h>

namespace cg {
template <typename Iter, typename Comp>
void sort_omp_async(Iter begin, Iter end, Comp comp) {
  size_t nthreads = omp_get_num_threads(), tid = omp_get_thread_num();
  size_t size = std::distance(begin, end);

  auto sort_beg = begin + (size * tid) / nthreads,
       sort_end = begin + (size * (tid + 1)) / nthreads;
  std::sort(sort_beg, sort_end, comp);
#pragma omp barrier

  for (size_t k = 1; k < nthreads; k *= 2) {
    auto beg_idx = (size * k * tid) / nthreads;
    auto mid_idx = (size * k * (tid + 1)) / nthreads;
    if (mid_idx < size) {
      auto end_idx = (size * k * (tid + 2)) / nthreads;
      end_idx = std::min(end_idx, size);

      auto merge_beg = begin + beg_idx, merge_mid = begin + mid_idx,
           merge_end = begin + end_idx;
      std::inplace_merge(merge_beg, merge_mid, merge_end, comp);
    }
#pragma omp barrier
  }
}

template <typename Iter>
void sort_omp_async(Iter begin, Iter end) {
  using T = typename std::iterator_traits<Iter>::value_type;
  sort_omp_async(begin, end, std::less<T>());
}
} // namespace cg