#pragma once
#include "../bit/view.h"
#include "../def/view.h"
#include "const_iterator.h"
#include "iterator.h"
#include "lane_ref.h"
#include "ref.h"
#include "sparse_ref.h"

namespace nitro::ind {
template <bool Indexed, typename T>
struct _view_impl;

template <typename T>
using view = typename _view_impl<is_indexed_v<T>, T>::type;

template <typename T>
struct _view_impl<false, T> {
  using type = def::view<T>;
};

// template <>
// struct _view_impl<false, bool> {
//   using type = bit::view;
// };

template <typename T>
struct _view_impl<true, T> {
  template <typename Idxes>
  struct _1;

  template <std::size_t... Idxes>
  struct _1<ind_seq<Idxes...>> {
    template <typename Types>
    struct _2;

    template <typename... Types>
    struct _2<type_list<Types...>> {
      class impl : public tuple<view<Types>...> {
      public:
        using tuple<view<Types>...>::tuple;

        impl() : tuple<view<Types>...>{view<Types>()...} {}

        ref<T> operator[](int idx) const {
          return ref<T>(this->template get<Idxes>()[idx]...);
        }

        template <typename Idx,
                  typename = std::enable_if_t<def::is_lane_like_v<Idx>>>
        auto operator[](Idx idx) const {
          return sparse_ref<T, Idx>(this->template get<Idxes>()[idx]...);
        }

        template <typename Idx>
        decltype(auto) at(Idx idx) const {
          return (*this)[idx];
        }

        template <std::size_t N, std::size_t W = def::opt_width_v>
        lane_ref<T, N, W> at_lane(int idx) const {
          return lane_ref<T, N, W>(
              this->template get<Idxes>().template at_lane<N, W>(idx)...);
        }

        int size() const {
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

        iterator<T> begin() const {
          return iterator<T>(this->template get<Idxes>().begin()...);
        }

        iterator<T> end() const {
          return iterator<T>(this->template get<Idxes>().end()...);
        }

        inline operator const_view<T>() const {
          return const_view<T>(
              static_cast<const_view<Types>>(this->template get<Idxes>())...);
        }

      private:
        using tuple<view<Types>...>::get;
      };
    };
  };

  using type = typename _1<idxes_t<T>>::template _2<subtypes_t<T>>::impl;
};
} // namespace nitro::ind