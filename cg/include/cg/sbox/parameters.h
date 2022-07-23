#pragma once
#include <cg/files/files.h>
#include <cg/types/amp.h>
#include <cg/utils/quantity.h>
#include <variant>

namespace cg::sbox {
struct parameters {
  quantity target_vel, accel_time, rest_period, avg_forces_over;

  struct attr_walls_t {
    std::string when, type;
    void load(ioxx::xyaml::node const &n);
  };
  attr_walls_t attr_walls;

  struct squeezing_t {
    bool perform;
    quantity target_density, vel_above_2V, accel_time;
    void load(ioxx::xyaml::node const &n);
  };
  squeezing_t squeezing;

  struct force_min_t {
    bool perform;
    quantity force_for_max_vel;
    void load(ioxx::xyaml::node const &n);
  };
  force_min_t force_min;

  struct oscillations_t {
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

    void load(ioxx::xyaml::node const &n);
  };
  oscillations_t oscillations;

  struct walls_t {
    std::string x_axis, y_axis, z_axis;
    quantity threshold;

    struct solid_wall_t {
      quantity depth;
      void load(ioxx::xyaml::node const &n);
    };
    solid_wall_t solid_wall;

    struct lj_wall_t {
      quantity depth, cycle_dur, breaking_dist_factor;
      std::optional<int> limit;
      void load(ioxx::xyaml::node const &n);
    };
    lj_wall_t lj_wall;

    struct harmonic_wall_t {
      quantity HH1;
      std::optional<int> limit;
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

  struct pulling_at_the_end_t {
    bool perform;
    void load(ioxx::xyaml::node const &n);
  };
  pulling_at_the_end_t pulling_at_the_end;

  void load(ioxx::xyaml::node const &n);
};
} // namespace cg::sbox