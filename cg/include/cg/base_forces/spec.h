#pragma once
#include "harmonic.h"
#include "lj.h"
#include "shifted_lj.h"
#include "sink_lj.h"
#include <cg/files/files.h>
#include <optional>
#include <variant>

namespace cg {
struct force_spec {
  std::string variant;
  std::optional<harmonic_specs> harmonic;
  std::optional<lj_specs> lj;
  std::optional<shifted_lj_specs> shifted_lj;
  std::optional<sink_lj_specs> sink_lj;

  void load(ioxx::xyaml::node const &node);
};

struct ss_force_spec {
  std::string variant;
  std::optional<ss_lj_specs> lj;
  std::optional<ss_shifted_lj_specs> shifted_lj;
  std::optional<ss_sink_lj_specs> sink_lj;

  void load(ioxx::xyaml::node const& node);
};
} // namespace cg