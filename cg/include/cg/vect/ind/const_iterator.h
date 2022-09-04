#pragma once
#include "../bit/const_iterator.h"
#include "../def/const_iterator.h"
#include "const_ref.h"

namespace nitro::ind {
template <bool Indexed, typename T>
struct _const_iterator_impl;

template <typename T>
using const_iterator = typename _const_iterator_impl<is_indexed_v<T>, T>::type;

template <typename T>
struct _const_iterator_impl<false, T> {
  using type = def::const_iterator<T>;
};

// template <>
// struct _const_iterator_impl<false, bool> {
//   using type = bit::const_iterator;
// };

template <typename T, typename Idxes>
class _ind_const_iterator;

template <typename T, std::size_t... Idxes>
class _ind_const_iterator<T, ind_seq<Idxes...>>
    : public tuple<
          const_iterator<std::tuple_element_t<Idxes, subtypes_t<T>>>...> {
public:
  template <size_t I>
  using const_iterator_ =
      const_iterator<std::tuple_element_t<I, subtypes_t<T>>>;
  using Base = tuple<const_iterator_<Idxes>...>;

  using Base::Base;

  explicit _ind_const_iterator(const_ref<T> const &ref_)
      : Base(const_iterator_<Idxes>(ref_.template get<Idxes>())...) {}

  const_ref<T> operator*() const {
    return const_ref<T>(*this->template get<Idxes>()...);
  }

  auto &operator++() {
    (..., ++this->template get<Idxes>());
    return *this;
  }

  auto &operator+=(int offset) {
    (..., advance<Idxes>(offset));
    return *this;
  }

  auto &operator--() {
    (..., --this->template get<Idxes>());
    return *this;
  }

  auto &operator-=(int offset) {
    (..., backtrack<Idxes>(offset));
    return *this;
  }

  auto operator++(int) {
    auto sav = *this;
    ++*this;
    return sav;
  }

  auto operator+(int offset) const {
    auto res = *this;
    res += offset;
    return res;
  }

  auto operator-(int offset) const {
    auto res = *this;
    res -= offset;
    return res;
  }

  std::ptrdiff_t operator-(_ind_const_iterator const &other) const {
    return this->template get<0>() - other.template get<0>();
  }

  bool operator<(_ind_const_iterator const &other) const {
    return this->template get<0>() < other.template get<0>();
  }

  bool operator<=(_ind_const_iterator const &other) const {
    return this->template get<0>() <= other.template get<0>();
  }

  bool operator>(_ind_const_iterator const &other) const {
    return this->template get<0>() > other.template get<0>();
  }

  bool operator>=(_ind_const_iterator const &other) const {
    return this->template get<0>() >= other.template get<0>();
  }

  bool operator==(_ind_const_iterator const &other) const {
    return this->template get<0>() == other.template get<0>();
  }

  bool operator!=(_ind_const_iterator const &other) const {
    return this->template get<0>() != other.template get<0>();
  }

private:
  using Base::get;

  template <std::size_t I>
  void advance(int offset) {
    this->template get<I>() += offset;
  }

  template <std::size_t I>
  void backtrack(int offset) {
    this->template get<I>() -= offset;
  }
};

template <typename T>
struct _const_iterator_impl<true, T> {
  using type = _ind_const_iterator<T, idxes_t<T>>;
};

} // namespace nitro::ind

namespace std {
template <typename T, typename Idxes>
struct iterator_traits<nitro::ind::_ind_const_iterator<T, Idxes>> {
  using value_type = T;
  using difference_type = std::ptrdiff_t;
  using reference = nitro::ind::const_ref<T>;
  using iterator_category = std::random_access_iterator_tag;
};
} // namespace std