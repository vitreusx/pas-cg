#pragma once
// #include "../bit/vector.h"
#include "../def/cuda_vector.h"
#include "allocator.h"
#include "const_iterator.h"
#include "const_view.h"
#include "cuda_allocator.h"
#include "iterator.h"
#include "vector.h"
#include "view.h"

namespace nitro::ind {
template <bool Indexed, typename T, typename Alloc>
struct _cuda_vector_impl;

template <typename T, typename Alloc = cuda_allocator<T>>
using cuda_vector = typename _cuda_vector_impl<is_indexed_v<T>, T, Alloc>::type;

template <typename T, typename Alloc>
struct _cuda_vector_impl<false, T, Alloc> {
  using type = def::cuda_vector<T, Alloc>;
};

// template <typename Alloc>
// struct _vector_impl<false, bool, Alloc> {
//   using type = bit::vector<Alloc>;
// };

template <typename T, typename Alloc>
struct _cuda_vector_impl<true, T, Alloc> {
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
      using _cuda_vector = cuda_vector<_subtype<I>, _alloc<I>>;

      using Base = tuple<_cuda_vector<Idxes>...>;

      class impl : public Base {
      public:
        impl() : Base{_cuda_vector<Idxes>()...} {}

        explicit impl(Alloc alloc)
            : Base{_cuda_vector<Idxes>(alloc.template get<Idxes>())...} {}

        template <typename CPUAlloc = allocator<T>>
        void pull_from(vector<T, CPUAlloc> const &v, cudaStream_t stream = 0) {
          (..., this->template get<Idxes>().pull_from(v.template get<Idxes>(),
                                                      stream));
        }

        template <typename CPUAlloc = allocator<T>>
        void push_to(vector<T, CPUAlloc> &v, cudaStream_t stream = 0) const {
          (..., this->template get<Idxes>().push_to(v.template get<Idxes>(),
                                                    stream));
        }

        void clear() {
          (..., this->template get<Idxes>().clear());
        }

        void reserve(int new_capacity) {
          (..., this->template get<Idxes>().reserve(new_capacity));
        }

        __host__ __device__ int size() const {
          return this->template get<0>().size();
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