#pragma once
#include "../def/array.h"
#include "const_iterator.h"
#include "const_view.h"
#include "iterator.h"
#include "view.h"

namespace nitro::ind {
template <bool Indexed, typename T, std::size_t Nx>
struct _array_impl;

template <typename T, std::size_t Nx>
using array = typename _array_impl<is_indexed_v<T>, T, Nx>::type;

template <typename T, std::size_t Nx>
struct _array_impl<false, T, Nx> {
  using type = def::array<T, Nx>;
};

// template <typename Alloc>
// struct _vector_impl<false, bool, Alloc> {
//   using type = bit::vector<Alloc>;
// };

template <typename T, std::size_t Nx>
struct _array_impl<true, T, Nx> {
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
      using _array = array<_subtype<I>, Nx>;

      using Base = tuple<_array<Idxes>...>;

      class impl : public Base {
      public:
        template <typename... Args>
        impl(Args &&...args)
            : Base{_ctor<Idxes, Args...>(std::forward<Args>(args)...)...} {} 

        int size() const {
          return this->template get<0>().size();
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

      private:
        using Base::get;

        template <std::size_t I, typename... Args>
        _array<I> _ctor(Args &&...args) {
          return _array<I>(args.template get<I>()...);
        }
      };
    };
  };

  using type = typename _1<idxes_t<T>>::template _2<subtypes_t<T>>::impl;
};
} // namespace nitro::ind