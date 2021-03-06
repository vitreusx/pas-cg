//#include <cg/afm/stats.h>
//#include <cg/utils/quantity.h>
// namespace cg::afm {
//
// void report_stats::report_to(out::report_data &report) const {
//  using namespace ioxx::xyaml;
//
//  if (report.report_files) {
//    csv<> vel_afm_csv;
//    vel_afm_csv.path = "vel-afm.csv";
//    vel_afm_csv.data.header = {"res idx", "force [f77unit*A/tau**2]"};
//    for (auto const &tip : tips->vel) {
//      auto &row = vel_afm_csv.data.emplace_back();
//      row["res idx"] = tip.res_idx();
//      auto force = eval_vel_forces->compute_force(tip);
//      row["force [f77unit*A/tau**2]"] =
//      quantity(force).in("f77unit*A/tau**2");
//    }
//
//    report.for_snap["const-velocity AFM"] << vel_afm_csv;
//
//    csv<> force_afm_csv;
//    force_afm_csv.path = "force-afm.csv";
//    force_afm_csv.data.header = {"res idx", "vel [A/tau]"};
//    for (auto const &tip : tips->force) {
//      auto &row = force_afm_csv.data.emplace_back();
//      row["res idx"] = tip.res_idx();
//      auto vel = norm(v[tip.res_idx()]);
//      row["vel [A/tau]"] = quantity(vel).in("A/tau");
//    }
//
//    report.for_snap["const-force AFM"] << force_afm_csv;
//
//    csv<> pulled_chains_csv;
//    pulled_chains_csv.path = "pulled-chains.csv";
//    pulled_chains_csv.data.header = {"chain idx", "end-to-end distance [A]"};
//    for (auto const &chain : tips->pulled_chains) {
//      auto &row = pulled_chains_csv.data.emplace_back();
//      row["chain idx"] = chain;
//      auto first = chain_first[chain], last = chain_last[chain];
//      auto dist = norm(r[last] - r[first]);
//      row["end-to-end distance [A]"] = quantity(dist).in("A");
//    }
//
//    report.for_snap["pulled-apart chains"] << pulled_chains_csv;
//  }
//
//  if (tips->pulled_chains.size() == 1) {
//    auto chain = tips->pulled_chains[0];
//    auto first = chain_first[chain], last = chain_last[chain];
//    auto dist = norm(r[last] - r[first]);
//    report.add_step_scalar("chain length [A]", quantity(dist).in("A"));
//  }
//}
//} // namespace cg::afm