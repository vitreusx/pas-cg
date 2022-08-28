#include "cg/utils/math.h"
#include <cg/simul/runtime.h>

namespace cg::simul {
task_slice::task_slice(cg::simul::sliceable_task const *task, int from, int to)
    : task{task}, from{from}, to{to} {}

void task_slice::perform() const {
  task->for_slice(from, to);
}

void sliceable_task::emplace_slices(std::vector<task_slice> &slices,
                                    int slice_size) const {
  auto total_size_ = total_size();
  for (int from = 0; from < total_size_;) {
    int to = min(from + slice_size, total_size_);
    slices.emplace_back(this, from, to);
    from = to;
  }
}

set_of_task_slices::set_of_task_slices(int num_threads)
    : num_threads{num_threads} {}

void set_of_task_slices::add_static_task(
    const cg::simul::sliceable_task &task) {
  static_range_tasks.push_back(&task);
  reset_all();
}

void set_of_task_slices::add_dynamic_task(
    const cg::simul::sliceable_task &task) {
  dynamic_range_tasks.push_back(&task);
  reset_dynamic();
}

void set_of_task_slices::reset_dynamic() {
  slices.erase(slices.begin() + num_fixed, slices.end());
  for (auto const *task : dynamic_range_tasks) {
    auto slice_size = (task->total_size() + num_threads - 1) / num_threads;
    task->emplace_slices(slices, slice_size);
  }
}

void set_of_task_slices::reset_all() {
  slices.clear();
  
  for (auto const *task : static_range_tasks) {
    auto slice_size = (task->total_size() + num_threads - 1) / num_threads;
    task->emplace_slices(slices, slice_size);
  }

  for (auto const *task : dynamic_range_tasks) {
    auto slice_size = (task->total_size() + num_threads - 1) / num_threads;
    task->emplace_slices(slices, slice_size);
  }

  num_fixed = (int)slices.size();
}

void set_of_task_slices::run_async(int thread_id) const {
  for (int task_idx = thread_id; task_idx < (int)slices.size();
       task_idx += num_threads)
    slices[task_idx].perform();
}
} // namespace cg::simul