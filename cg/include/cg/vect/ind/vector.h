#pragma once
#include "../bit/vector.h"
#include "../def/vector.h"
#include "allocator.h"
#include "const_iterator.h"
#include "const_view.h"
#include "iterator.h"
#include "view.h"

namespace nitro::ind {
template <bool Indexed, typename T, typename Alloc>
struct _vector_impl;

template <typename T, typename Alloc = allocator<T>>
using vector = typename _vector_impl<is_indexed_v<T>, T, Alloc>::type;

template <typename T, typename Alloc>
struct _vector_impl<false, T, Alloc> {
  using type = def::vector<T, Alloc>;
};

// template <typename Alloc>
// struct _vector_impl<false, bool, Alloc> {
//   using type = bit::vector<Alloc>;
// };

template <typename T, typename Alloc>
struct _vector_impl<true, T, Alloc> {
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
      using _vector = vector<_subtype<I>, _alloc<I>>;

      using Base = tuple<_vector<Idxes>...>;

      class impl : public Base {
      public:
        impl() : Base{_vector<Idxes>()...} {}

        explicit impl(Alloc alloc)
            : Base{_vector<Idxes>(alloc.template get<Idxes>())...} {}

        explicit impl(int n, T const &init = T(), Alloc alloc = Alloc())
            : Base{_vector<Idxes>(n, init.template get<Idxes>(),
                                  alloc.template get<Idxes>())...} {}

        template <typename E, typename = std::enable_if_t<is_expr_for<T, E>>>
        explicit impl(int n, E const &e, Alloc alloc = Alloc())
            : Base{_vector<Idxes>(n, e.template get<Idxes>(),
                                  alloc.template get<Idxes>())...} {}

        int size() const {
          return this->template get<0>().size();
        }

        int capacity() const {
          return this->template get<0>().capacity();
        }

        ref<T> operator[](int idx) {
          return ref<T>(this->template get<Idxes>()[idx]...);
        }

        const_ref<T> operator[](int idx) const {
          return const_ref<T>(this->template get<Idxes>()[idx]...);
        }

        ref<T> at(int idx) {
          return (*this)[idx];
        }

        const_ref<T> at(int idx) const {
          return (*this)[idx];
        }

        template <std::size_t N, std::size_t W>
        lane_ref<T, N, W> at_lane(int idx) {
          return lane_ref<T, N, W>(
              this->template get<Idxes>().template at_lane<N, W>(idx)...);
        }

        template <std::size_t N, std::size_t W>
        lane_const_ref<T, N, W> at_lane(int idx) const {
          return lane_const_ref<T, N, W>(
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

      private:
        using Base::get;
      };
    };
  };

  using type = typename _1<idxes_t<T>>::template _2<subtypes_t<T>>::impl;
};
} // namespace nitro::ind