#pragma once
#include "../def/dual_vector.h"
#include "cuda_allocator.h"
#include "vector.h"

namespace nitro::ind {
template <bool Indexed, typename T, typename Alloc, typename CudaAlloc>
struct _dual_vector_impl;

template <typename T, typename Alloc = allocator<T>,
          typename CudaAlloc = cuda_allocator<T>>
using dual_vector =
    typename _dual_vector_impl<is_indexed_v<T>, T, Alloc, CudaAlloc>::type;

template <typename T, typename Alloc, typename CudaAlloc>
struct _dual_vector_impl<false, T, Alloc, CudaAlloc> {
  using type = def::dual_vector<T, Alloc, CudaAlloc>;
};

// template <typename Alloc>
// struct _vector_impl<false, bool, Alloc> {
//   using type = bit::vector<Alloc>;
// };

template <typename T, typename Alloc, typename CudaAlloc>
struct _dual_vector_impl<true, T, Alloc, CudaAlloc> {
  template <typename Idxes>
  struct _1;

  template <std::size_t... Idxes>
  struct _1<ind_seq<Idxes...>> {
    template <typename Types>
    struct _2;

    template <typename... Types>
    struct _2<type_list<Types...>> {
      using _Types = type_list<Types...>;

      template <std::size_t I>
      using _subtype = std::tuple_element_t<I, _Types>;

      template <std::size_t I>
      using _alloc = std::tuple_element_t<I, subtypes_t<Alloc>>;

      template <std::size_t I>
      using _cuda_alloc = std::tuple_element_t<I, subtypes_t<CudaAlloc>>;

      template <std::size_t I>
      using _dual_vector = dual_vector<_subtype<I>, _alloc<I>, _cuda_alloc<I>>;

      using Base = tuple<_dual_vector<Idxes>...>;

      class impl : public Base {
      public:
        impl() : Base{_dual_vector<Idxes>()...} {}

        explicit impl(Alloc alloc, CudaAlloc cuda_alloc)
            : Base{_dual_vector<Idxes>(alloc.template get<Idxes>(),
                                       cuda_alloc.template get<Idxes>())...} {}

        explicit impl(int n, T const &init = T(), Alloc alloc = Alloc(),
                      CudaAlloc cuda_alloc = CudaAlloc())
            : Base{_dual_vector<Idxes>(n, init.template get<Idxes>(),
                                       alloc.template get<Idxes>())...} {}

        template <typename E, typename = std::enable_if_t<is_expr_for<T, E>>>
        explicit impl(int n, E const &e, Alloc alloc = Alloc(),
                      CudaAlloc cuda_alloc = CudaAlloc())
            : Base{_dual_vector<Idxes>(n, e.template get<Idxes>(),
                                       alloc.template get<Idxes>(),
                                       cuda_alloc.template get<Idxes>())...} {}

        int size() const {
          return this->template get<0>().size();
        }

        int capacity() const {
          return this->template get<0>().capacity();
        }

        ref<T> operator[](int idx) {
          return ref<T>(this->template get<Idxes>()[idx]...);
        }

        template <typename Idx,
                  typename = std::enable_if_t<def::is_lane_like_v<Idx>>>
        auto operator[](Idx idx) {
          return sparse_ref<T, Idx>(this->template get<Idxes>()[idx]...);
        }

        template <typename Idx, typename Mask,
                  typename = std::enable_if_t<def::is_lane_like_v<Idx> &&
                                              def::is_lane_like_v<Mask>>>
        auto operator[](std::pair<Idx, Mask> idx_mask) {
          return masked_ref<T, Idx, Mask>{
              this->template get<Idxes>()[idx_mask]...};
        }

        const_ref<T> operator[](int idx) const {
          return const_ref<T>(this->template get<Idxes>()[idx]...);
        }

        template <typename Idx,
                  typename = std::enable_if_t<def::is_lane_like_v<Idx>>>
        auto operator[](Idx idx) const {
          return sparse_const_ref<T, Idx>(this->template get<Idxes>()[idx]...);
        }

        template <typename Idx, typename Mask,
                  typename = std::enable_if_t<def::is_lane_like_v<Idx> &&
                                              def::is_lane_like_v<Mask>>>
        auto operator[](std::pair<Idx, Mask> idx_mask) const {
          return masked_const_ref<T, Idx, Mask>{
              this->template get<Idxes>()[idx_mask]...};
        }

        template <typename Idx>
        decltype(auto) at(Idx const &idx) {
          return (*this)[idx];
        }

        template <typename Idx>
        decltype(auto) at(Idx const &idx) const {
          return (*this)[idx];
        }

        template <std::size_t N, std::size_t W = def::opt_width_v>
        lane_ref<T, N, W> at_lane(int idx) {
          return lane_ref<T, N, W>(
              this->template get<Idxes>().template at_lane<N, W>(idx)...);
        }

        template <std::size_t N, std::size_t W = def::opt_width_v>
        lane<T, N, W> at_lane(int idx) const {
          return lane<T, N, W>(
              this->template get<Idxes>().template at_lane<N, W>(idx)...);
        }

        template <std::size_t N>
        int num_lanes() const {
          return size() / N;
        }

        template <std::size_t N>
        int final_idx() const {
          return N * num_lanes<N>();
        }

        void clear() {
          (..., this->template get<Idxes>().clear());
        }

        void reserve(int new_capacity) {
          (..., this->template get<Idxes>().reserve(new_capacity));
        }

        void resize(int new_size, T const &init = T()) {
          (..., this->template get<Idxes>().resize(new_size,
                                                   init.template get<Idxes>()));
        }

        template <typename E, typename = std::enable_if_t<is_expr_for<T, E>>>
        void resize(int new_size, E const &e) {
          (..., this->template get<Idxes>().resize(new_size,
                                                   e.template get<Idxes>()));
        }

        void shrink(int new_size) {
          (..., this->template get<Idxes>().shrink(new_size));
        }

        template <typename E, typename = std::enable_if_t<is_expr_for<T, E>>>
        void push_back(E const &e) {
          (..., this->template get<Idxes>().push_back(e.template get<Idxes>()));
        }

        template <typename... Args>
        ref<T> emplace_back(Args &&...args) {
          push_back(T(std::forward<Args>(args)...));
          return (*this)[size() - 1];
        }

        iterator<T> begin() {
          return iterator<T>(this->template get<Idxes>().begin()...);
        }

        const_iterator<T> begin() const {
          return const_iterator<T>(this->template get<Idxes>().begin()...);
        }

        iterator<T> end() {
          return iterator<T>(this->template get<Idxes>().end()...);
        }

        const_iterator<T> end() const {
          return const_iterator<T>(this->template get<Idxes>().end()...);
        }

        operator view<T>() {
          return view<T>(static_cast<view<_subtype<Idxes>>>(
              this->template get<Idxes>())...);
        }

        operator const_view<T>() const {
          return const_view<T>(static_cast<const_view<_subtype<Idxes>>>(
              this->template get<Idxes>())...);
        }

      public:
        void push_to_gpu() {
          (..., this->template get<Idxes>().push_to_gpu());
        }

        void pull_from_gpu() {
          (..., this->template get<Idxes>().pull_from_gpu());
        }

        view<T> gpu_view() {
          return view<T>(this->template get<Idxes>().gpu_view()...);
        }

        const_view<T> gpu_view() const {
          return view<T>(this->template get<Idxes>().gpu_view()...);
        }

      private:
        using Base::get;
      };
    };
  };

  using type = typename _1<idxes_t<T>>::template _2<subtypes_t<T>>::impl;
};
} // namespace nitro::ind