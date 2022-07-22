#include <cg/sbox/parameters.h>
#include <iostream>

namespace cg::sbox {
void parameters::squeezing_t::load(const ioxx::xyaml::node &n) {
  n["perform"] >> perform;
  n["target density"] >> target_density;
  n["velocity above 2V"] >> vel_above_2V;
  n["velocity below 2V"] >> vel_below_2V;
  n["acceleration time"] >> accel_time;
}

void parameters::force_min_t::load(const ioxx::xyaml::node &n) {
  n["perform"] >> perform;
  n["max velocity"] >> max_velocity;
  n["force for max velocity"] >> force_for_max_vel;
}

void parameters::oscilations_t::amplitude_t::load(const ioxx::xyaml::node &n) {
  n["variant"] >> variant;
  n["absolute value"] >> abs_value;
  n["relative value"] >> rel_value;
}

void parameters::oscilations_t::load(const ioxx::xyaml::node &n) {
  n["perform"] >> perform;
  n["type"] >> type;
  n["num of cycles"] >> num_of_cycles;
  n["amplitude"] >> amplitude;
  n["angular frequency"] >> angular_freq;

  auto an = n["attractive walls"];
  an["when"] >> when_attr;
  an["type"] >> attr_wall_type;
}

void parameters::walls_t::solid_wall_t::load(const ioxx::xyaml::node &n) {
  n["depth"] >> depth;
}

void parameters::walls_t::lj_wall_t::load(const ioxx::xyaml::node &n) {
  n["depth"] >> depth;
  n["cycle duration"] >> cycle_dur;
  n["breaking dist factor"] >> breaking_dist_factor;
}

void parameters::walls_t::harmonic_wall_t::load(const ioxx::xyaml::node &n) {
  n["HH1"] >> HH1;
}

void parameters::walls_t::load(const ioxx::xyaml::node &n) {
  if (n["all axes"]) {
    auto val = n["all axes"].as<std::string>();
    x_axis = y_axis = z_axis = val;
  }

  if (n["x axis"])
    n["x axis"] >> x_axis;
  if (n["y axis"])
    n["y axis"] >> y_axis;
  if (n["z axis"])
    n["z axis"] >> z_axis;

  n["threshold"] >> threshold;
  n["solid wall params"] >> solid_wall;
  n["lj wall params"] >> lj_wall;
  n["harmonic wall params"] >> harmonic_wall;
}

void parameters::init_size_t::sufficient_t::load(const ioxx::xyaml::node &n) {
  n["pad (x bond)"] >> pad;
  n["maximum density"] >> max_density;
  n["cubic"] >> cubic;
}

void parameters::init_size_t::load(const ioxx::xyaml::node &n) {
  n["type"] >> type;
  n["sufficient box params"] >> sufficient;
}

void parameters::pulling_at_the_end_t::load(const ioxx::xyaml::node &n) {
  n["perform"] >> perform;
  n["velocity"] >> velocity;
}

void parameters::load(const ioxx::xyaml::node &n) {
  n["initial size"] >> init_size;
  n["squeezing"] >> squeezing;
  n["finding force minimum"] >> force_min;
  n["oscilations"] >> oscilations;
  n["walls"] >> walls;
  n["rest period"] >> rest_period;
  n["pulling at the end"] >> pulling_at_the_end;
}
} // namespace cg::sbox