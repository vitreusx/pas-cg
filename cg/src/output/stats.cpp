#include "output/stats.h"
#include "utils/math.h"
using namespace cg::out;

void add_stats::report_to(report_state &report) const {
  auto stats_node = report.current["stats"];
  stats_node["t"] = *t;
  stats_node["V"] = *V;

  real K = 0;
  for (int idx = 0; idx < v.size(); ++idx)
    K += (real)0.5 * mass[(uint8_t)atype[idx]] * ipow<2>(norm(v[idx]));

  stats_node["K"] = K;
  stats_node["E"] = K + *V;
}