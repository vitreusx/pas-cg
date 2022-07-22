#pragma once
#include <cg/files/files.h>
#include <cg/types/amp.h>
#include <cg/utils/quantity.h>
#include <variant>

namespace cg::sbox {
struct parameters {
  struct squeezing_t {
    bool perform;
    quantity target_density;
    std::string wall_type;
    void load(ioxx::xyaml::node const &n);
  };
  squeezing_t squeezing;

  struct movement_t {
    bool perform;
    std::string type;
    int num_of_cycles;

    struct amplitude_t {
      std::string variant;
      quantity abs_value, rel_value;
      void load(ioxx::xyaml::node const &n);
    };
    amplitude_t amplitude;

    quantity angular_freq;
    std::string when_attr, attr_wall_type;

    void load(ioxx::xyaml::node const &n);
  };
  movement_t movement;

  struct walls_t {
    std::string x_axis, y_axis, z_axis;
    quantity threshold;

    struct solid_wall_t {
      quantity depth, cutoff;
      void load(ioxx::xyaml::node const &n);
    };
    solid_wall_t solid_wall;

    struct lj_wall_t {
      quantity depth, cycle_dur, breaking_dist_factor;
      void load(ioxx::xyaml::node const &n);
    };
    lj_wall_t lj_wall;

    struct harmonic_wall_t {
      quantity HH1, breaking_dist_factor;
      void load(ioxx::xyaml::node const &n);
    };
    harmonic_wall_t harmonic_wall;

    void load(ioxx::xyaml::node const &n);
  };
  walls_t walls;

  struct init_size_t {
    std::string type;

    struct sufficient_t {
      quantity pad;
      std::optional<quantity> max_density;
      bool cubic;
      void load(ioxx::xyaml::node const &n);
    };
    sufficient_t sufficient;

    void load(ioxx::xyaml::node const &n);
  };
  init_size_t init_size;

  void load(ioxx::xyaml::node const &n);
};
} // namespace cg::sbox