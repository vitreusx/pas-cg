#pragma once

namespace nitro::def {
template <typename T, typename Repr>
inline constexpr bool is_representable_as_v = sizeof(T) == sizeof(Repr) &&
                                              alignof(T) == alignof(Repr);

template <typename T>
struct type_identity {
  using type = T;
};

template <typename T, typename... Reprs>
struct _repr_impl;

template <typename T>
struct _repr_impl<T> {};

template <typename T, typename Head, typename... Tail>
struct _repr_impl<T, Head, Tail...>
    : std::conditional_t<is_representable_as_v<T, Head>, type_identity<Head>,
                       _repr_impl<T, Tail...>> {};

template <typename T>
using _repr = _repr_impl<T, uint8_t, uint16_t, uint32_t, uint64_t>;

template <typename T>
using repr_t = typename _repr<T>::type;
} // namespace nitro::def