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
  vector() : data{nullptr}, _size{0}, _capacity{0}, alloc{Alloc()} {}

  explicit vector(Alloc alloc)
      : data{nullptr}, _size{0}, _capacity{0}, alloc{alloc} {};

  explicit vector(int n, T const &init = T(), Alloc alloc = Alloc()) {
    this->alloc = alloc;
    data = nullptr;
    _size = _capacity = n;

    if (n > 0) {
      data = this->alloc.allocate(n);
      std::uninitialized_fill_n(data, n, init);
    }
  }

  ~vector() {
    destroy();
  }

  vector(vector const &other) {
    alloc = other.alloc;
    data = nullptr;
    _size = _capacity = other._size;

    if (other._size > 0) {
      data = alloc.allocate(other._size);
      std::uninitialized_copy_n(other.data, other._size, data);
    }
  }

  vector &operator=(vector const &other) {
    if (this != &other) {
      destroy();
      alloc = other.alloc;
      if (other._size > 0) {
        data = alloc.allocate(other._size);
        std::uninitialized_copy_n(other.data, other._size, data);
      }
      _size = _capacity = other._size;
    }
    return *this;
  }

  vector(vector &&other) noexcept {
    data = other.data;
    _size = other._size;
    _capacity = other._capacity;

    other.data = nullptr;
    other._size = other._capacity = 0;
  }

  vector &operator=(vector &&other) noexcept {
    if (this != &other) {
      destroy();
      data = other.data;
      _size = other._size;
      _capacity = other._capacity;
      other.data = nullptr;
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
    return data[idx];
  }

  template <typename Idxes, typename = std::enable_if_t<is_lane_like_v<Idxes>>>
  auto operator[](Idxes idxes) {
    return sparse_ref<T, Idxes>(data, idxes);
  }

  T const &operator[](int idx) const {
    return data[idx];
  }

  template <typename Idxes, typename = std::enable_if_t<is_lane_like_v<Idxes>>>
  auto operator[](Idxes idxes) const {
    using Data = lane<T, lane_size_v<Idxes>, lane_width_v<Idxes>>;
    return gather<Data>(data, idxes);
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
    return lane_ref<T, N, W>(data + N * idx);
  }

  template <std::size_t N, std::size_t W = opt_width_v>
  lane<T, N, W> at_lane(int idx) const {
    return construct<lane<T, N, W>>(data + N * idx);
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
    if (data && _size > 0)
      std::destroy_n(data, _size);

    _size = 0;
  }

  void reserve(int new_capacity) {
    if (new_capacity <= _capacity)
      return;

    auto new_data = alloc.allocate(new_capacity);
    if (data) {
      if (_size > 0) {
        std::uninitialized_move_n(data, _size, new_data);
        std::destroy_n(data, _size);
      }
      alloc.deallocate(data, _capacity);
    }
    data = new_data;
    _capacity = new_capacity;
  }

  void resize(int new_size, T const &init = T()) {
    if (new_size < _size) {
      std::destroy_n(data + new_size, _size - new_size);
    } else if (new_size > _size) {
      reserve(new_size);
      std::uninitialized_fill_n(data + _size, new_size - _size, init);
    }
    _size = new_size;
  }

  void shrink(int new_size) {
    if (new_size < _size) {
      std::destroy_n(data + new_size, _size - new_size);
      _size = new_size;
    }
  }

  void push_back(T const &value) {
    if (_size + 1 > _capacity)
      reserve(2 * _size + 8);

    std::uninitialized_copy_n(&value, 1, data + _size);
    ++_size;
  }

  template <typename... Args>
  T &emplace_back(Args &&...args) {
    if (_size >= _capacity)
      reserve(2 * (_size + 1) + 8);

    ::new (data + _size) T(std::forward<Args>(args)...);
    return data[_size++];
  }

  auto begin() {
    return iterator<T>(data);
  }

  auto begin() const {
    return const_iterator<T>(data);
  }

  auto end() {
    return iterator<T>(data + size());
  }

  auto end() const {
    return const_iterator<T>(data + size());
  }

  operator view<T>() {
    return view<T>(data, size());
  }

  operator const_view<T>() const {
    return const_view<T>(data, size());
  }

private:
  void destroy() {
    if (data) {
      if (_size > 0)
        std::destroy_n(data, _size);
      alloc.deallocate(data, _capacity);
      data = nullptr;
    }
    _size = _capacity = 0;
  }

  T *data;
  int _size, _capacity;
  Alloc alloc;
};
} // namespace nitro::def