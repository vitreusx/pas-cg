#pragma once
#include <vector>

namespace cg::simul {
class sliceable_task;

class task_slice {
public:
  explicit task_slice(sliceable_task const *task, int from, int to);
  void perform() const;

private:
  sliceable_task const *task;
  int from, to;
};

class sliceable_task {
public:
  virtual void for_slice(int from, int to) const = 0;
  void emplace_slices(std::vector<task_slice> &slices, int slice_size) const;
  virtual int total_size() const = 0;
};

class set_of_task_slices {
public:
  set_of_task_slices() = default;
  set_of_task_slices(int num_threads);

  void add_static_task(sliceable_task const &task);
  void add_dynamic_task(sliceable_task const &task);

  void reset_all();
  void reset_dynamic();

  void run_async(int thread_id) const;

private:
  int num_fixed = 0, num_threads = 1;
  std::vector<task_slice> slices;
  std::vector<sliceable_task const *> static_range_tasks, dynamic_range_tasks;
};
} // namespace cg::simul