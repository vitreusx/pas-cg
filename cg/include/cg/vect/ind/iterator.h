#pragma once
#include "../bit/iterator.h"
#include "../def/iterator.h"
#include "ref.h"

namespace nitro::ind {
template <bool Indexed, typename T>
struct _iterator_impl;

template <typename T>
using iterator = typename _iterator_impl<is_indexed_v<T>, T>::type;

template <typename T>
struct _iterator_impl<false, T> {
  using type = def::iterator<T>;
};

// template <>
// struct _iterator_impl<false, bool> {
//   using type = bit::iterator;
// };

template <typename T, typename Idxes>
class _ind_iterator;

template <typename T, std::size_t... Idxes>
class _ind_iterator<T, ind_seq<Idxes...>>
    : public tuple<iterator<std::tuple_element_t<Idxes, subtypes_t<T>>>...> {
public:
  template <size_t I>
  using iterator_ = iterator<std::tuple_element_t<I, subtypes_t<T>>>;
  using Base = tuple<iterator_<Idxes>...>;

  using Base::Base;

  explicit _ind_iterator(ref<T> const &ref_)
      : Base(iterator_<Idxes>(ref_.template get<Idxes>())...) {}

  ref<T> operator*() const {
    return ref<T>(*this->template get<Idxes>()...);
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

  std::ptrdiff_t operator-(_ind_iterator const &other) const {
    return this->template get<0>() - other.template get<0>();
  }

  bool operator<(_ind_iterator const &other) const {
    return this->template get<0>() < other.template get<0>();
  }

  bool operator<=(_ind_iterator const &other) const {
    return this->template get<0>() <= other.template get<0>();
  }

  bool operator>(_ind_iterator const &other) const {
    return this->template get<0>() > other.template get<0>();
  }

  bool operator>=(_ind_iterator const &other) const {
    return this->template get<0>() >= other.template get<0>();
  }

  bool operator==(_ind_iterator const &other) const {
    return this->template get<0>() == other.template get<0>();
  }

  bool operator!=(_ind_iterator const &other) const {
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
struct _iterator_impl<true, T> {
  using type = _ind_iterator<T, idxes_t<T>>;
};

} // namespace nitro::ind

namespace std {
template <typename T, typename Idxes>
struct iterator_traits<nitro::ind::_ind_iterator<T, Idxes>> {
  using value_type = T;
  using difference_type = std::ptrdiff_t;
  using reference = nitro::ind::ref<T>;
  using iterator_category = std::random_access_iterator_tag;
};
} // namespace std