#include "output/stats.h"
#include "utils/ioxx_interop.h"
#include "utils/math.h"
#include <ioxx/ioxx.h>
using namespace cg::out;

struct chain_row {
  int chain_idx, first, last;
  cg::real end_to_end_dist;

  void connect(ioxx::row_proxy &row) {
    row["idx"] << chain_idx;
    row["first"] << first;
    row["last"] << last;
    row["end to end dist[A]"] << end_to_end_dist;
  }
};

void add_stats::report_to(report_state &report) const {
  auto stats_node = report.current["stats"];
  stats_node["t"] = *t;
  stats_node["V"] = *V;

  real K = 0;
  for (int idx = 0; idx < v.size(); ++idx)
    K += (real)0.5 * mass[(uint8_t)atype[idx]] * ipow<2>(norm(v[idx]));

  stats_node["K"] = K;
  stats_node["E"] = K + *V;

  ioxx::xyaml::csv<chain_row> chains_csv;
  chains_csv.path = "chains.csv";
  chains_csv.data.header = {"idx", "first", "last", "end to end dist[A]"};

  int num_chains = chain_first.size();
  for (int idx = 0; idx < num_chains; ++idx) {
    chain_row row{};
    row.chain_idx = idx;
    row.first = chain_first[idx];
    row.last = chain_last[idx];
    auto dist = norm(r[row.last] - r[row.first]);
    row.end_to_end_dist = quantity(dist).in("A");
    chains_csv.data.rows.push_back(row);
  }
  report.current["chains"] = chains_csv;
}