#include "pbar/render.h"
#include "utils/text.h"
#include <iostream>
namespace cg::pbar {

static std::string timedelta(uint64_t ms) {
  int ms_part = ms % 1000;
  uint64_t rem_s = (ms - ms_part) / 1000;
  int s_part = rem_s % 60;
  uint64_t rem_min = (rem_s - s_part) / 60;
  int min_part = rem_min % 60;
  uint64_t rem_hr = (rem_min - min_part) / 60;
  int hr_part = rem_hr % 24;
  int day_part = (rem_hr - hr_part) / 24;

  return cg::format("%d:%02d:%02d:%02d.%04d", day_part, hr_part, min_part,
                    s_part, ms_part);
}

void render::operator()() const {
  using namespace std;
  using namespace std::chrono;

  if (*is_first) {
    *start_clock = *last_clock = high_resolution_clock::now();
  }

  auto now = high_resolution_clock::now();
  auto since_last = duration_cast<nanoseconds>(now - *last_clock);
  auto since_last_s = (real)since_last.count() / 1e9;

  if (*is_first || since_last_s >= period_s) {
    auto since_start = duration_cast<milliseconds>(now - *start_clock);
    auto since_start_ms = since_start.count();

    auto progress = *t / total_time;

    cout << "\r"
         << "[";
    auto pos = (int)floor(width * progress);
    for (int i = 0; i < width; ++i) {
      if (i < pos)
        cout << "=";
      else if (i == pos)
        cout << ">";
      else
        cout << " ";
    }
    cout << "] " << *t << " / " << total_time;
    cout << " V = " << *V;
    cout << " t = " << timedelta(since_start_ms);

    cout.flush();

    *last_clock = now;
    *is_first = false;
  }
}}