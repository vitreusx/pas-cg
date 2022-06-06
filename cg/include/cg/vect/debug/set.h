#pragma once
#include "../common/set_node.h"
#include "vector.h"

namespace nitro::debug {
template <typename T> using set_node = nitro::set_node<T>;

template <typename T, typename Alloc = std::allocator<set_node<T>>>
class set : public vector<set_node<T>, Alloc> {
public:
  using Base = vector<set_node<T>, Alloc>;

  set() : Base() {}

  explicit set(Alloc const &alloc) : Base(alloc) {}

  explicit set(int n, T const &init = T(), Alloc const &alloc = Alloc())
      : Base(n, set_node<T>(init, false), alloc) {}

  template <typename U>
  explicit set(int n, U const &init, Alloc const &alloc = Alloc())
      : Base(n, set_node<T>(init, false), alloc) {}

  template <typename U> void insert(U const &value) {
    push_back(set_node<T>(value, false));
  }

  template <typename... Args> set_node<T> &emplace(Args &&...args) {
    return vector<set_node<T>, Alloc>::emplace_back(
        T(std::forward<Args>(args)...), false);
  }

  void remove(int idx) {
    (*this)[idx].remove();
  }

  int num_items() const {
    int res = 0;
    for (int idx = 0; idx < this->size(); ++idx) {
      if ((*this)[idx].has_item())
        ++res;
    }
    return res;
  }

  void compress() {
    int slot_idx = 0, num_items = 0;
    for (int cur_idx = 0; cur_idx < this->size(); ++cur_idx) {
      if ((*this)[cur_idx].has_item()) {
        (*this)[slot_idx++].item() = (*this)[cur_idx].item();
        ++num_items;
      }
    }
    this->shrink(num_items);
  }
};
} // namespace nitro::debug