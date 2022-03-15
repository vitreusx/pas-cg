#pragma once
#include <cg/types/amp.h>

namespace cg {
template <typename E> struct sync_data_expr : public nitro::ind_expr<E> {
  EXPR_BODY(back, side_all, side_polar, side_hydrophobic)

  template <typename F>
  inline auto &operator+=(sync_data_expr<F> const &other) {
    back() += other.back();
    side_all() += other.side_all();
    side_polar() += other.side_polar();
    side_hydrophobic() += other.side_hydrophobic();
    return *this;
  }

  template <typename F>
  inline auto &operator-=(sync_data_expr<F> const &other) {
    back() -= other.back();
    side_all() -= other.side_all();
    side_polar() -= other.side_polar();
    side_hydrophobic() -= other.side_hydrophobic();
    return *this;
  }

  inline bool is_valid() const {
    return (back() >= 0) && (side_all() >= 0) && (side_polar() >= 0) &&
           (side_hydrophobic() >= 0);
  }
};

template <typename E1, typename E2>
class sync_sum_expr : public sync_data_expr<sync_sum_expr<E1, E2>> {
public:
  sync_sum_expr(E1 const &e1, E2 const &e2) : e1{e1}, e2{e2} {};

  int8_t back() const { return e1.back() + e2.back(); }

  int8_t side_all() const { return e1.side_all() + e2.side_all(); }

  int8_t side_polar() const { return e1.side_polar() + e2.side_polar(); }

  int8_t side_hydrophobic() const {
    return e1.side_hydrophobic() + e2.side_hydrophobic();
  }

private:
  E1 e1;
  E2 e2;
};

template <typename E1, typename E2>
inline auto operator+(sync_data_expr<E1> const &e1,
                      sync_data_expr<E2> const &e2) {
  return sync_sum_expr<E1, E2>(static_cast<E1 const &>(e1),
                               static_cast<E2 const &>(e2));
}

template <typename E1, typename E2>
class sync_diff_expr : public sync_data_expr<sync_diff_expr<E1, E2>> {
public:
  sync_diff_expr(E1 const &e1, E2 const &e2) : e1{e1}, e2{e2} {};

  int8_t back() const { return e1.back() - e2.back(); }

  int8_t side_all() const { return e1.side_all() - e2.side_all(); }

  int8_t side_polar() const { return e1.side_polar() - e2.side_polar(); }

  int8_t side_hydrophobic() const {
    return e1.side_hydrophobic() - e2.side_hydrophobic();
  }

private:
  E1 e1;
  E2 e2;
};

template <typename E1, typename E2>
inline auto operator-(sync_data_expr<E1> const &e1,
                      sync_data_expr<E2> const &e2) {
  return sync_diff_expr<E1, E2>(static_cast<E1 const &>(e1),
                                static_cast<E2 const &>(e2));
}

template <typename E> struct sync_data_auto_expr : public sync_data_expr<E> {
  AUTO_EXPR_BODY(back, side_all, side_polar, side_hydrophobic)
};

using sync_data_base = nitro::tuple_wrapper<int8_t, int8_t, int8_t, int8_t>;

class sync_data : public sync_data_auto_expr<sync_data>, public sync_data_base {
public:
  using Base = sync_data_base;
  using Base::Base;
  using Base::get;

  sync_data() : Base(0, 0, 0, 0){};
};
} // namespace cg

namespace nitro {
template <> struct is_indexed_impl<cg::sync_data> : public std::true_type {};

template <typename E> struct expr_impl<E, cg::sync_data> {
  using type = cg::sync_data_expr<E>;
};

template <typename E> struct auto_expr_impl<E, cg::sync_data> {
  using type = cg::sync_data_auto_expr<E>;
};
}; // namespace nitro