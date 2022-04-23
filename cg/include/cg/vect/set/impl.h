#pragma once
#include "../vector.h"
#include "node.h"

namespace nitro {
template <typename T, typename Alloc = allocator<set_node<T>>,
          typename Idx = int>
class set : public vector<set_node<T>, Alloc, Idx> {
public:
  using Base = vector<set_node<T>, Alloc, Idx>;

  set() : Base() {}

  explicit set(Alloc const &alloc) : Base(alloc) {}

  explicit set(Idx n, T const &init = T(), Alloc const &alloc = Alloc())
      : Base(n, set_node<T>(init, false), alloc) {}

  template <typename U>
  explicit set(Idx n, U const &init, Alloc const &alloc = Alloc())
      : Base(n, set_node<T>(init, false), alloc) {}

  template <typename U> void insert(U const &value) {
    push_back(set_node<T>(value, false));
  }

  template <typename... Args> at_expr<set_node<T>> emplace(Args &&...args) {
    return emplace_back(T(std::forward<Args>(args)...), false);
  }

  void remove(Idx idx) { (*this)[idx].remove(); }

  Idx num_items() const {
    Idx res = 0;
    for (Idx idx = 0; idx < this->size(); ++idx) {
      if (this->at(idx).has_item())
        ++res;
    }
    return res;
  }

  void compress() {
    Idx slot_idx = 0, num_items = 0;
    for (Idx cur_idx = 0; cur_idx < this->size(); ++cur_idx) {
      if (this->at(cur_idx).has_item()) {
        this->at(slot_idx++).item() = this->at(cur_idx).item();
        ++num_items;
      }
    }
    this->shrink(num_items);
  }

private:
  using Base::emplace_back;
  using Base::push_back;
  using Base::reserve;
  using Base::resize;
};
} // namespace nitro