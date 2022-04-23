#pragma once

namespace nitro {
template <typename E, typename T> struct expr_impl;

template <typename E, typename T> using expr = typename expr_impl<E, T>::type;

template <typename E, typename T> struct auto_expr_impl;

template <typename E, typename T>
using auto_expr = typename auto_expr_impl<E, T>::type;

}
