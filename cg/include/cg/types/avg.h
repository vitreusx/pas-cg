#pragma once
#include <optional>

namespace cg {
template <typename ValueT, typename TimeT> class moving_avg {
public:
  moving_avg() = default;
  explicit moving_avg(TimeT const &window) : window{window} {}

  std::optional<ValueT> const *operator->() const {
    return &avg;
  }

  void add(TimeT const &time, ValueT const &sample) {
    if (num_samples == 0) {
      init_time = time;
      sum = sample;
      num_samples = 1;
    } else {
      sum += sample;
      ++num_samples;
    }

    if (init_time + window <= time) {
      avg = sum / num_samples;
      num_samples = 0;
    }
  }

private:
  std::optional<ValueT> avg;
  ValueT sum;
  TimeT init_time = 0, window = 0;
  int num_samples = 0;
};
} // namespace cg