#pragma once
#include "../lane.h"
#include "decl.h"

namespace nitro {
template <typename T, size_t N> struct def_lane_at {
public:
  explicit def_lane_at(T &ref) : ref{ref} {}
  def_lane_at(def_lane_at &other) = default;
  def_lane_at(def_lane_at &&other) noexcept = default;

  operator lane<T, N>() {
    lane<T, N> res;
    res.load(&ref);
    return res;
  }

  def_lane_at &operator=(def_lane_at const &other) {
    assign(other);
    return *this;
  }

  template <typename E> def_lane_at &operator=(E const &e) {
    assign(e);
    return *this;
  }

  def_lane_at &operator=(def_lane_at &&other) noexcept {
    assign(other);
    return *this;
  }

  template <typename E> def_lane_at &operator=(E &&e) noexcept {
    assign(e);
    return *this;
  }

private:
  template <typename E> void assign(E const &e) noexcept {
    lane<T, N> x = e;
    x.store(&ref);
  }

  T &ref;
};
} // namespace nitro