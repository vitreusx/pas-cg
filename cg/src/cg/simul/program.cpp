#include <cg/simul/program.h>
#include <iostream>
#include <optional>
#include <string>
namespace cg::simul {

static std::vector<std::string> split(std::string const &str, char sep) {
  std::vector<std::string> parts = {""};
  for (auto ch : str) {
    if (ch == sep) {
      parts.emplace_back("");
    } else {
      parts.back().push_back(ch);
    }
  }
  return parts;
}

struct program_args {
  std::string prog;
  bool help = false;
  std::optional<std::string> error;
  std::vector<std::filesystem::path> param_paths;
  std::optional<std::filesystem::path> ckpt_path;
  std::vector<std::string> setters;

  enum class state {
    PROG_NAME,
    FIRST_ARG,
    OPTION,
    PARAM_PATH_VALUE,
    HELP,
    INVALID,
    CKPT_PATH_VALUE,
    SET_VALUE,
  };

  void print_help();
  program_args(int argc, char **argv);
};

void program_args::print_help() {
  std::cout << "PAS CG Coarse Grained Model" << '\n';
  std::cout << "Documentation: https://vitreusx.github.io/pas-cg" << '\n';
  std::cout << "Usage: " << prog << " [options...]" << '\n';
  std::cout << "  -h, --help                                       Print this "
               "help message."
            << '\n';
  std::cout << "  -p [path], --param-file [path]                   Load "
               "parameter file."
            << '\n';
  std::cout << "  -c [path], --ckpt-file [path]                    Load "
               "simulation from checkpoint."
            << '\n';
  std::cout << "  -s \"key_1.key_2.(...).key_n=value\", --set ...  Override "
               "input file setting."
            << '\n';
}

program_args::program_args(int argc, char **argv) {
  state parser_state = state::PROG_NAME;
  for (int idx = 0; idx < argc; ++idx) {
    auto arg = std::string(argv[idx]);
    switch (parser_state) {
    case state::PROG_NAME:
      prog = arg;
      parser_state = state::FIRST_ARG;
      break;
    case state::FIRST_ARG:
    case state::OPTION:
      if (arg == "-h" || arg == "--help") {
        help = true;
        parser_state = state::HELP;
      } else if (arg == "-p" || arg == "--param-file") {
        parser_state = state::PARAM_PATH_VALUE;
      } else if (arg == "-c" || arg == "--ckpt-file") {
        parser_state = state::CKPT_PATH_VALUE;
      } else if (arg == "-s" || arg == "--set") {
        parser_state = state::SET_VALUE;
      } else {
        std::stringstream err_ss;
        err_ss << "Invalid option \"" << arg << "\"";
        error = err_ss.str();
        parser_state = state::INVALID;
      }
      break;
    case state::PARAM_PATH_VALUE:
      param_paths.emplace_back(arg);
      parser_state = state::OPTION;
      break;
    case state::CKPT_PATH_VALUE:
      if (!ckpt_path.has_value()) {
        ckpt_path = arg;
        parser_state = state::OPTION;
      } else {
        error = "Can provide only one checkpoint file.";
        parser_state = state::INVALID;
      }
      break;
    case state::SET_VALUE:
      setters.push_back(arg);
      parser_state = state::OPTION;
      break;
    case state::HELP:
    case state::INVALID:
      break;
    }
  }

  if (parser_state == state::FIRST_ARG)
    help = true;
}

void program::main(int argc, char **argv) {
  auto args = program_args(argc, argv);

  if (args.error.has_value() || args.help) {
    if (args.error.has_value()) {
      std::cout << "Error: " << args.error.value() << '\n';
    }

    args.print_help();
    auto exit_code = args.help ? EXIT_SUCCESS : EXIT_FAILURE;
    exit(exit_code);
  } else if (args.ckpt_path.has_value()) {
    run_from_checkpoint(args.ckpt_path.value());
  } else {
    parameters params;
    using namespace ioxx::xyaml;
    auto params_yml = defaults_yml();
    for (auto const &path : args.param_paths) {
      auto slice_yml = node::import(path);
      params_yml.merge(slice_yml);
    }

    for (auto const &setter : args.setters) {
      auto p1 = split(setter, '=');
      auto key = p1[0], value = p1[1];

      std::vector<ioxx::xyaml::node> descent = {params_yml};
      for (auto const &key_ : split(key, '.'))
        descent.push_back(descent.back()[key_]);
      descent.back() = value;
    }

    params_yml >> params;

    using prog_mode = gen::parameters::prog_mode;
    switch (params.gen.mode) {
    case prog_mode::perform_simulation:
      perform_simulation(params, params_yml);
      break;
    case prog_mode::check_determinism:
      check_determinism(params, params_yml);
      break;
    }
  }
}

void program::perform_simulation(parameters const &params,
                                 ioxx::xyaml::node const &raw_params) {
  omp_set_num_threads((int)params.gen.num_of_threads);

  auto st = state(params, raw_params);
  auto team = thread_team(st);

#pragma omp parallel default(none) shared(team)
  {
    auto &thr = team.fork();
    thr.main();
  }
}

void program::check_determinism(parameters const &params,
                                ioxx::xyaml::node const &raw_params) {
  omp_set_num_threads((int)params.gen.num_of_threads);

  auto st1 = state(params, raw_params), st2 = state(params, raw_params);
  auto team1 = thread_team(st1), team2 = thread_team(st2);

#pragma omp parallel default(none) shared(st1, team1, st2, team2)
  {
    auto &thr1 = team1.fork(), &thr2 = team2.fork();

    do {
      thr1.step();
#pragma omp barrier
      thr2.step();
#pragma omp barrier

#pragma omp master
      st1.verify_equal(st2);
#pragma omp barrier
    } while (st1.cur_phase != phase::SIMUL_END &&
             st2.cur_phase != phase::SIMUL_END);
  }
}

void program::run_from_checkpoint(const std::filesystem::path &ckpt_path) {
  std::ifstream ckpt_is(ckpt_path, std::ios::binary);
  state st;
  ckpt_is >> st;

  auto const &params = st.params;
  omp_set_num_threads((int)params.gen.num_of_threads);

  auto team = thread_team(st);

#pragma omp parallel default(none) shared(team)
  {
    auto &thr = team.fork();
    thr.main();
  }
}
} // namespace cg::simul