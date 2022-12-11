#pragma once
#include <tuple>
#include <vector>

namespace cg::simul {
class task {
public:
  virtual ~task() = default;
  virtual void run() const = 0;
};

template <typename Iterable>
class iter_task_mixin : public task {
public:
  void run() const override {
    auto n = static_cast<Iterable const &>(*this).size();
    for (int idx = 0; idx < n; ++idx)
      static_cast<Iterable const &>(*this).iter(idx);
  }

  void omp_async() const {
    auto n = static_cast<Iterable const &>(*this).size();
#pragma omp for schedule(static) nowait
    for (int idx = 0; idx < n; ++idx)
      static_cast<Iterable const &>(*this).iter(idx);
  }
};

template <typename Iterable>
class vect_iter_task_mixin : public task {
public:
  void run() const override {
    static constexpr auto V = Iterable::elems_per_vect;
    auto n = static_cast<Iterable const &>(*this).size(), vect_n = n / V;

    for (int vect_idx = 0; vect_idx < vect_n; ++vect_idx)
      static_cast<Iterable const &>(*this).vect_iter(vect_idx);

    for (int idx = vect_n * V; idx < n; ++idx)
      static_cast<Iterable const &>(*this).iter(idx);
  }

  void omp_async() const {
    static constexpr auto V = Iterable::elems_per_vect;
    auto n = static_cast<Iterable const &>(*this).size(), vect_n = n / V;

#pragma omp for schedule(static) nowait
    for (int vect_idx = 0; vect_idx < vect_n; ++vect_idx)
      static_cast<Iterable const &>(*this).vect_iter(vect_idx);

#pragma omp for schedule(static) nowait
    for (int idx = vect_n * V; idx < n; ++idx)
      static_cast<Iterable const &>(*this).iter(idx);
  }
};

class set_of_tasks : public task {
public:
  std::vector<task const *> subtasks;

public:
  void run() const override;
  void omp_async() const;
};

class divisible {
public:
  virtual ~divisible() = default;
  virtual void reset() = 0;
  virtual void divide(std::vector<task const *> &subtasks, int size_hint) = 0;
};

class set_of_divisibles : public task {
public:
  std::vector<std::tuple<divisible *, bool const *, int>> divs;

  void update();
  void run() const override;
  void omp_async() const;

private:
  set_of_tasks tasks;
};

template <typename Iterable>
class iter_divisible_mixin : public divisible,
                             public iter_task_mixin<Iterable> {
public:
  void reset() override {
    slices.clear();
  }

  void divide(std::vector<task const *> &subtasks, int size_hint) override {
    auto n = static_cast<Iterable const &>(*this).size();
    auto this_ = static_cast<Iterable const *>(this);
    for (int from = 0; from < n; from += size_hint) {
      auto to = std::min(from + size_hint, n);
      slices.emplace_back(this_, from, to);
    }

    for (auto& slice_: slices)
      subtasks.push_back(&slice_);
  }

private:
  class slice : public task {
  public:
    explicit slice(Iterable const *iterable, int from, int to)
        : iterable{iterable}, from{from}, to{to} {}

    void run() const override {
      for (int idx = from; idx < to; ++idx)
        iterable->iter(idx);
    }

  private:
    Iterable const *iterable;
    int from, to;
  };

  std::vector<slice> slices;
};

template <typename Iterable>
class vect_iter_divisible_mixin : public divisible,
                                  public vect_iter_task_mixin<Iterable> {
public:
  void reset() override {
    slices.clear();
    vect_slices.clear();
  }

  void divide(std::vector<task const *> &subtasks, int size_hint) override {
    static constexpr auto V = Iterable::elems_per_vect;
    auto n = static_cast<Iterable const &>(*this).size(), vect_n = n / V;
    int vect_size_hint = size_hint / V;
    auto this_ = static_cast<Iterable const *>(this);

    for (int vect_from = 0; vect_from < vect_n; vect_from += vect_size_hint) {
      int vect_to = std::min(vect_from + vect_size_hint, vect_n);
      vect_slices.emplace_back(this_, vect_from, vect_to);
    }

    for (auto& slice_: vect_slices)
      subtasks.push_back(&slice_);

    for (int from = vect_n * V; from < n; from += size_hint) {
      int to = std::min(from + size_hint, n);
      slices.emplace_back(this_, from, to);
    }

    for (auto& slice_: slices)
      subtasks.push_back(&slice_);
  }

private:
  class slice : public task {
  public:
    explicit slice(Iterable const *iterable, int from, int to)
        : iterable{iterable}, from{from}, to{to} {}

    void run() const override {
      for (int idx = from; idx < to; ++idx)
        iterable->iter(idx);
    }

  private:
    Iterable const *iterable;
    int from, to;
  };

  class vect_slice : public task {
  public:
    explicit vect_slice(Iterable const *iterable, int from, int to)
        : iterable{iterable}, from{from}, to{to} {}

    void run() const override {
      for (int idx = from; idx < to; ++idx)
        iterable->vect_iter(idx);
    }

  private:
    Iterable const *iterable;
    int from, to;
  };

  std::vector<slice> slices;
  std::vector<vect_slice> vect_slices;
};

} // namespace cg::simul