#pragma once
// #include "../bit/const_view.h"
#include "../def/const_view.h"
#include "const_iterator.h"
#include "const_ref.h"
#include "lane.h"
#include "masked_const_ref.h"
#include "sparse_const_ref.h"

namespace nitro::ind {
template <bool Indexed, typename T>
struct _const_view_impl;

template <typename T>
using const_view = typename _const_view_impl<is_indexed_v<T>, T>::type;

template <typename T>
struct _const_view_impl<false, T> {
  using type = def::const_view<T>;
};

// template <>
// struct _const_view_impl<false, bool> {
//   using type = bit::const_view;
// };

template <typename T>
struct _const_view_impl<true, T> {
  template <typename Idxes>
  struct _1;

  template <std::size_t... Idxes>
  struct _1<ind_seq<Idxes...>> {
    template <typename Types>
    struct _2;

    template <typename... Types>
    struct _2<type_list<Types...>> {
      class impl : public tuple<const_view<Types>...> {
      public:
        using tuple<const_view<Types>...>::tuple;

        impl() : tuple<const_view<Types>...>{const_view<Types>()...} {}

        __host__ __device__ const_ref<T> operator[](int idx) const {
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
        __host__ __device__ decltype(auto) at(Idx idx) const {
          return (*this)[idx];
        }

        template <std::size_t N, std::size_t W>
        lane<T, N, W> at_lane(int idx) const {
          return lane<T, N, W>(
              this->template get<Idxes>().template at_lane<N, W>(idx)...);
        }

        __host__ __device__ int size() const {
          return this->template get<0>().size();
        }

        template <std::size_t N>
        int num_lanes() const {
          return this->template get<0>().template num_lanes<N>();
        }

        template <std::size_t N>
        int final_idx() const {
          return N * num_lanes<N>();
        }

        const_iterator<T> begin() const {
          return const_iterator<T>(this->template get<Idxes>().begin()...);
        }

        const_iterator<T> end() const {
          return const_iterator<T>(this->template get<Idxes>().end()...);
        }

      private:
        using tuple<const_view<Types>...>::get;
      };
    };
  };

  using type = typename _1<idxes_t<T>>::template _2<subtypes_t<T>>::impl;
};
} // namespace nitro::ind