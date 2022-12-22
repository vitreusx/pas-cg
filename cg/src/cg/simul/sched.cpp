#include <cg/simul/sched.h>

namespace cg::simul {
void set_of_tasks::run() const {
  for (auto const *subtask : subtasks)
    subtask->run();
}

void set_of_tasks::omp_async() const {
#pragma omp for schedule(static, 1) nowait
  // #pragma omp for schedule(dynamic, 1) nowait
  for (auto const *subtask : subtasks)
    subtask->run();
}

void set_of_divisibles::update() {
  tasks.subtasks.clear();
  for (auto const &[div_task, is_enabled, size_hint] : divs) {
    if (*is_enabled) {
      div_task->reset();
      div_task->divide(tasks.subtasks, size_hint);
    }
  }
}

void set_of_divisibles::run() const {
  tasks.run();
}

void set_of_divisibles::omp_async() const {
  tasks.omp_async();
}
} // namespace cg::simul