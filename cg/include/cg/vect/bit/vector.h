#pragma once
#include "allocator.h"
#include "const_iterator.h"
#include "const_view.h"
#include "iterator.h"
#include "view.h"
#include <vector>

namespace nitro::bit {
template <typename Alloc = allocator>
class vector {
public:
  vector()
      : data{nullptr}, _size{0}, _byte_size{0},
        _byte_capacity{0}, alloc{Alloc()} {}

  explicit vector(Alloc alloc)
      : data{nullptr}, _size{0}, _byte_size{0},
        _byte_capacity{0}, alloc{alloc} {};

  explicit vector(int n, bool const &init = bool(), Alloc alloc = Alloc()) {
    this->alloc = alloc;
    data = nullptr;
    _size = n;
    _byte_size = _byte_capacity = req_bytes(n);

    if (_byte_size > 0) {
      data = this->alloc.allocate(_size);
      auto byte_init = init ? ~(byte)0 : (byte)0;
      std::uninitialized_fill_n(data, _byte_size, byte_init);
    }
  }

  ~vector() {
    destroy();
  }

  vector(vector const &other) {
    alloc = other.alloc;
    data = nullptr;
    _size = other._size;
    _byte_size = _byte_capacity = other._byte_size;

    if (other._size > 0) {
      data = alloc.allocate(other._size);
      std::uninitialized_copy_n(other.data, other._byte_size, data);
    }
  }

  vector &operator=(vector const &other) {
    if (&*this != &other) {
      destroy();
      alloc = other.alloc;
      if (other._size > 0) {
        data = alloc.allocate(other._size);
        std::uninitialized_copy_n(other.data, other._byte_size, data);
      }
      _size = other._size;
      _byte_size = _byte_capacity = other._byte_size;
    }
    return *this;
  }

  vector(vector &&other) noexcept {
    data = other.data;
    _size = other._size;
    _byte_size = other._byte_size;
    _byte_capacity = other._byte_capacity;

    other.data = nullptr;
    other._size = other._byte_size = other._byte_capacity = 0;
  }

  vector &operator=(vector &&other) noexcept {
    if (&*this != &other) {
      destroy();
      data = other.data;
      _size = other._size;
      _byte_size = other._byte_size;
      _byte_capacity = other._byte_capacity;

      other.data = nullptr;
      other._size = other._byte_size = other._byte_capacity = 0;
    }
    return *this;
  }

  int size() const {
    return _size;
  }

  int capacity() const {
    return 8 * _byte_capacity;
  }

  ref operator[](int idx) {
    return ref(data + byte_offset(idx), byte_mask(idx));
  }

  const_ref operator[](int idx) const {
    return const_ref(data + byte_offset(idx), byte_mask(idx));
  }

  ref at(int idx) {
    return (*this)[idx];
  }

  const_ref at(int idx) const {
    return (*this)[idx];
  }

  template <std::size_t N, std::size_t W>
  lane_ref<N, W> at_lane(int idx) {
    return lane_ref<N, W>(data + byte_offset(N * idx));
  }

  template <std::size_t N, std::size_t W>
  lane_const_ref<N, W> at_lane(int idx) const {
    return lane_const_ref<N, W>(data + byte_offset(N * idx));
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
      std::destroy_n(data, _byte_size);

    _size = 0;
  }

  void reserve(int new_capacity) {
    auto new_byte_capacity = req_bytes(new_capacity);
    if (new_byte_capacity <= _byte_capacity)
      return;

    auto new_data = alloc.allocate(new_capacity);
    if (data) {
      if (_size > 0) {
        std::uninitialized_move_n(data, _byte_size, new_data);
        std::destroy_n(data, _byte_size);
      }
      alloc.deallocate(data, _size);
    }
    data = new_data;
    _byte_capacity = new_byte_capacity;
  }

  void resize(int new_size, bool const &init = bool()) {
    auto new_byte_size = req_bytes(new_size);
    if (new_byte_size < _byte_size) {
      std::destroy_n(data + new_size, _size - new_size);
    } else if (new_byte_size > _byte_size) {
      reserve(new_byte_size);
      auto byte_init = init ? ~(byte)0 : (byte)0;
      std::uninitialized_fill_n(data + _byte_size, new_byte_size - _byte_size,
                                byte_init);
    }
    _size = new_size;
  }

  void shrink(int new_size) {
    if (new_size < _size) {
      auto new_byte_size = req_bytes(new_size);
      if (new_byte_size < _byte_size)
        std::destroy_n(data + new_size, _size - new_size);
      _size = new_size;
    }
  }

  void push_back(bool const &value) {
    auto new_size = _size + 1, new_byte_size = (int)req_bytes(new_size);
    if (new_byte_size > _byte_capacity)
      reserve(2 * _size + 8);

    (*this)[_size] = value;
    _size = new_size;
    _byte_size = new_byte_size;
  }

  template <typename... Args>
  ref emplace_back(Args &&...args) {
    auto new_size = _size + 1, new_byte_size = (int)req_bytes(new_size);
    if (new_byte_size > _byte_capacity)
      reserve(2 * _size + 8);

    bool val;
    ::new (&val) bool(std::forward<Args>(args)...);
    auto ref = (*this)[_size];
    ref = val;

    _size = new_size;
    _byte_size = new_byte_size;
    return ref;
  }

  auto begin() {
    return iterator(data);
  }

  auto begin() const {
    return const_iterator(data);
  }

  auto end() {
    return iterator(data, size());
  }

  auto end() const {
    return const_iterator(data, size());
  }

  operator view() {
    return view(data, size());
  }

  operator const_view() const {
    return const_view(data, size());
  }

private:
  void destroy() {
    if (data) {
      if (_size > 0)
        std::destroy_n(data, _byte_size);
      alloc.deallocate(data, _size);
      data = nullptr;
    }
    _size = _byte_capacity = 0;
  }

  byte *data;
  int _size, _byte_size, _byte_capacity;
  Alloc alloc;
};
} // namespace nitro::bit