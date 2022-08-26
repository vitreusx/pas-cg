#include "cg/utils/math.h"
#include <cg/simul/runtime.h>

namespace cg::simul {
task_slice::task_slice(cg::simul::sliceable_task const *task, int from, int to)
    : task{task}, from{from}, to{to} {}

void task_slice::perform() const {
  task->for_slice(from, to);
}

void sliceable_task::emplace_slices(std::vector<task_slice> &slices) const {
  auto total_size_ = total_size(), slice_size_ = slice_size();
  for (int from = 0; from < total_size_;) {
    int to = min(from + slice_size_, total_size_);
    slices.emplace_back(this, from, to);
    from = to;
  }
}

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
  for (auto const *task : dynamic_range_tasks)
    task->emplace_slices(slices);
}

void set_of_task_slices::reset_all() {
  slices.clear();
  for (auto const *task : static_range_tasks)
    task->emplace_slices(slices);
  num_fixed = (int)slices.size();
  for (auto const *task : dynamic_range_tasks)
    task->emplace_slices(slices);
}

void set_of_task_slices::run_async(int thread_id, int num_threads) const {
  for (int task_idx = thread_id; task_idx < (int)slices.size();
       task_idx += num_threads)
    slices[task_idx].perform();
}
} // namespace cg::simul