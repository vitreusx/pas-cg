#pragma once

namespace nitro::aos {
template <bool IsIndexed, typename T>
struct allocator_impl;

template <typename T>
using allocator = typename allocator_impl<is_indexed_v<T>, T>::type;

template <typename T>
struct allocator_impl<false, T> {
  using type = std::allocator<T>;
};

template <typename Types>
class soa_allocator;

template <typename... Types>
class soa_allocator<type_list<Types...>> : public tuple<allocator<Types>...> {
public:
  soa_allocator() : tuple<allocator<Types>...>{allocator<Types>{}...} {
  }
};

template <typename T>
struct allocator_impl<true, T> {
  using type = soa_allocator<decltype(T::Types())>;
};
}