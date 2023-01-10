#pragma once
#include "allocator.h"
#include "const_iterator.h"
#include "const_view.h"
#include "iterator.h"
#include "view.h"

namespace nitro::def {
template <typename T, typename Alloc = allocator<T>>
class vector {
public:
  vector() : data_{nullptr}, _size{0}, _capacity{0}, alloc{Alloc()} {}

  explicit vector(Alloc alloc)
      : data_{nullptr}, _size{0}, _capacity{0}, alloc{alloc} {};

  explicit vector(int n, T const &init = T(), Alloc alloc = Alloc()) {
    this->alloc = alloc;
    data_ = nullptr;
    _size = _capacity = n;

    if (n > 0) {
      data_ = this->alloc.allocate(n);
      std::uninitialized_fill_n(data_, n, init);
    }
  }

  ~vector() {
    destroy();
  }

  vector(vector const &other) {
    alloc = other.alloc;
    data_ = nullptr;
    _size = _capacity = other._size;

    if (other._size > 0) {
      data_ = alloc.allocate(other._size);
      std::uninitialized_copy_n(other.data_, other._size, data_);
    }
  }

  vector &operator=(vector const &other) {
    if (this != &other) {
      destroy();
      alloc = other.alloc;
      if (other._size > 0) {
        data_ = alloc.allocate(other._size);
        std::uninitialized_copy_n(other.data_, other._size, data_);
      }
      _size = _capacity = other._size;
    }
    return *this;
  }

  vector(vector &&other) noexcept {
    data_ = other.data_;
    _size = other._size;
    _capacity = other._capacity;

    other.data_ = nullptr;
    other._size = other._capacity = 0;
  }

  vector &operator=(vector &&other) noexcept {
    if (this != &other) {
      destroy();
      data_ = other.data_;
      _size = other._size;
      _capacity = other._capacity;
      other.data_ = nullptr;
      other._size = other._capacity = 0;
    }
    return *this;
  }

  int size() const {
    return _size;
  }

  int capacity() const {
    return _capacity;
  }

  T &operator[](int idx) {
    return data_[idx];
  }

  template <typename Idxes, typename = std::enable_if_t<is_lane_like_v<Idxes>>>
  auto operator[](Idxes idxes) {
    return sparse_ref<T, Idxes>(data_, idxes);
  }

  template <typename Idxes, typename Mask,
            typename =
                std::enable_if_t<is_lane_like_v<Idxes> && is_lane_like_v<Mask>>>
  auto operator[](std::pair<Idxes, Mask> idxes_mask) {
    auto const &[idxes, mask] = idxes_mask;
    return masked_ref(data_, idxes, mask);
  }

  T const &operator[](int idx) const {
    return data_[idx];
  }

  template <typename Idxes, typename = std::enable_if_t<is_lane_like_v<Idxes>>>
  auto operator[](Idxes idxes) const {
    using Data = lane<T, lane_size_v<Idxes>, lane_width_v<Idxes>>;
    return gather<Data>(data_, idxes);
  }

  template <typename Idxes, typename Mask,
            typename =
                std::enable_if_t<is_lane_like_v<Idxes> && is_lane_like_v<Mask>>>
  auto operator[](std::pair<Idxes, Mask> idxes_mask) const {
    using Data = lane<T, lane_size_v<Idxes>, lane_width_v<Idxes>>;
    auto const &[idxes, mask] = idxes_mask;
    return masked_gather<Data>(data_, idxes, mask);
  }

  template <typename Idx>
  decltype(auto) at(Idx idx) {
    return (*this)[idx];
  }

  template <typename Idx>
  decltype(auto) at(Idx idx) const {
    return (*this)[idx];
  }

  template <std::size_t N, std::size_t W = opt_width_v>
  lane_ref<T, N, W> at_lane(int idx) {
    return lane_ref<T, N, W>(data_ + N * idx);
  }

  template <std::size_t N, std::size_t W = opt_width_v>
  lane<T, N, W> at_lane(int idx) const {
    return construct<lane<T, N, W>>(data_ + N * idx);
  }

  template <std::size_t N>
  int num_lanes() const {
    return size() / N;
  }

  template <std::size_t N>
  int final_idx() const {
    return N * num_lanes<N>();
  }

  void clear() {
    if (data_ && _size > 0)
      std::destroy_n(data_, _size);

    _size = 0;
  }

  void reserve(int new_capacity) {
    if (new_capacity <= _capacity)
      return;

    auto new_data = alloc.allocate(new_capacity);
    if (data_) {
      if (_size > 0) {
        std::uninitialized_move_n(data_, _size, new_data);
        std::destroy_n(data_, _size);
      }
      alloc.deallocate(data_, _capacity);
    }
    data_ = new_data;
    _capacity = new_capacity;
  }

  void resize(int new_size, T const &init = T()) {
    if (new_size < _size) {
      std::destroy_n(data_ + new_size, _size - new_size);
    } else if (new_size > _size) {
      reserve(new_size);
      std::uninitialized_fill_n(data_ + _size, new_size - _size, init);
    }
    _size = new_size;
  }

  void shrink(int new_size) {
    if (new_size < _size) {
      std::destroy_n(data_ + new_size, _size - new_size);
      _size = new_size;
    }
  }

  void push_back(T const &value) {
    if (_size + 1 > _capacity)
      reserve(2 * _size + 8);

    std::uninitialized_copy_n(&value, 1, data_ + _size);
    ++_size;
  }

  template <typename... Args>
  T &emplace_back(Args &&...args) {
    if (_size >= _capacity)
      reserve(2 * (_size + 1) + 8);

    ::new (data_ + _size) T(std::forward<Args>(args)...);
    return data_[_size++];
  }

  auto begin() {
    return iterator<T>(data_);
  }

  auto begin() const {
    return const_iterator<T>(data_);
  }

  auto end() {
    return iterator<T>(data_ + size());
  }

  auto end() const {
    return const_iterator<T>(data_ + size());
  }

  operator view<T>() {
    return view<T>(data_, size());
  }

  operator const_view<T>() const {
    return const_view<T>(data_, size());
  }

  T *data() {
    return data_;
  }

  T const *data() const {
    return data_;
  }

private:
  void destroy() {
    if (data_) {
      if (_size > 0)
        std::destroy_n(data_, _size);
      alloc.deallocate(data_, _capacity);
      data_ = nullptr;
    }
    _size = _capacity = 0;
  }

  T *data_;
  int _size, _capacity;
  Alloc alloc;
};
} // namespace nitro::def