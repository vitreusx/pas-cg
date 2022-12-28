#pragma once
#include "../cuda_interop.h"
#include "ind_seq.h"
#include "type_list.h"
#include "type_traits.h"
#include <utility>

namespace nitro::ind {
template <typename... Types>
class tuple;

template <>
class tuple<> {
public:
  static auto Idxes() {
    return std::make_index_sequence<0>{};
  }

  static auto Subtypes() {
    return type_list<>{};
  }
};

template <typename Head, typename... Tail>
class tuple<Head, Tail...> {
public:
  __host__ __device__ explicit tuple(Head &&head, Tail &&...tail)
      : head{std::forward<Head>(head)}, tail(std::forward<Tail>(tail)...) {}

  template <typename U = Head,
            typename = std::enable_if_t<std::conjunction_v<
                std::is_default_constructible<U>,
                std::is_default_constructible<tuple<Tail...>>>>>
  __host__ __device__ tuple() : head{}, tail{} {}

  template <typename E, typename = std::enable_if_t<is_indexed_v<E>>>
  __host__ __device__ auto &operator=(E const &e) {
    assign<E>(e, e.Idxes());
    return *this;
  }

  __host__ __device__ tuple(tuple const &) = default;
  __host__ __device__ tuple(tuple &&) = default;

  __host__ __device__ auto &operator=(tuple const &other) {
    assign<tuple>(other, other.Idxes());
    return *this;
  }

  template <std::size_t I>
  __host__ __device__ decltype(auto) get() {
    if constexpr (I == 0)
      return (head);
    else
      return tail.template get<I - 1>();
    __builtin_unreachable();
  }

  template <std::size_t I>
  __host__ __device__ decltype(auto) get() const {
    if constexpr (I == 0)
      return (head);
    else
      return tail.template get<I - 1>();
    __builtin_unreachable();
  }

  static auto Idxes() {
    return std::make_index_sequence<1 + sizeof...(Tail)>{};
  }

  static auto Subtypes() {
    return type_list<Head, Tail...>{};
  }

private:
  Head head;
  tuple<Tail...> tail;

  template <typename E, std::size_t I>
  __host__ __device__ void assign1(E const &e) {
    this->get<I>() = e.template get<I>();
  }

  template <typename E, std::size_t... _Idxes>
  __host__ __device__ void assign(E const &e, ind_seq<_Idxes...>) {
    (..., assign1<E, _Idxes>(e));
  }
};
} // namespace nitro::ind