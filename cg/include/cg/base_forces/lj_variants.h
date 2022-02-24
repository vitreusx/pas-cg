#pragma once
#include "lj.h"
#include "sink_lj.h"
#include <cg/amino/amino_acid.h>
#include <cg/utils/hash.h>
#include <ioxx/xyaml.h>
#include <unordered_map>

namespace cg {
struct lj_variants {
  lj bb, bs, sb;
  std::unordered_map<std::pair<amino_acid, amino_acid>, sink_lj> ss;

  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg