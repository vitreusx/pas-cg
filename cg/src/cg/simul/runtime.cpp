#include "cg/utils/math.h"
#include <cg/simul/runtime.h>

namespace cg::simul {
task_slice::task_slice(cg::simul::sliceable_task const *task, int from, int to)
    : task{task}, from{from}, to{to} {}

void task_slice::perform() const {
  task->for_slice(from, to);
}

set_of_task_slices::set_of_task_slices(int num_threads)
    : num_threads{num_threads} {}

void set_of_task_slices::add_static_task(const cg::simul::sliceable_task &task,
                                         const bool &enabled) {
  static_range_tasks.emplace_back(&task, &enabled);
}

void set_of_task_slices::add_dynamic_task(const cg::simul::sliceable_task &task,
                                          const bool &enabled) {
  dynamic_range_tasks.emplace_back(&task, &enabled);
}

void set_of_task_slices::emplace_from(
    const cg::simul::set_of_task_slices::task_list &tasks) {
  for (auto const &[task, enabled] : tasks) {
    if (!*enabled)
      continue;

    //    auto slice_size = (task->total_size() + num_threads - 1) /
    //    num_threads;
    auto slice_size = 256;
    emplace_slices(task, slice_size);
  }
}

void set_of_task_slices::reset_dynamic() {
  slices.erase(slices.begin() + num_fixed, slices.end());
  emplace_from(dynamic_range_tasks);
}

void set_of_task_slices::reset_all() {
  slices.clear();
  emplace_from(static_range_tasks);
  num_fixed = (int)slices.size();
  emplace_from(dynamic_range_tasks);
}

void set_of_task_slices::run_async(int thread_id) const {
  for (int task_idx = thread_id; task_idx < (int)slices.size();
       task_idx += num_threads)
    slices[task_idx].perform();
}

void set_of_task_slices::emplace_slices(sliceable_task const *task,
                                        int slice_size) {
  auto total_size_ = task->total_size();
  for (int from = 0; from < total_size_;) {
    int to = min(from + slice_size, total_size_);
    slices.emplace_back(task, from, to);
    from = to;
  }
}
} // namespace cg::simul