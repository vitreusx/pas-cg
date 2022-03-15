#include "simul/program.h"
#include <iostream>
#include <optional>
#include <string>
namespace cg::simul {

struct program_args {
  std::string prog;
  bool help = false;
  std::optional<std::string> error;
  std::vector<std::filesystem::path> param_paths;
  std::optional<std::filesystem::path> chkpt_path;

  enum class state {
    PROG_NAME,
    OPTION,
    PARAM_PATH_VALUE,
    CHKPT_PATH_VALUE,
    HELP,
    INVALID
  };

  void print_help();
  program_args(int argc, char **argv);
};

void program_args::print_help() {
  if (error.has_value()) {
    std::cout << "Error: " << error.value() << '\n';
  }

  std::cout << "Usage: " << prog << " [options...]" << '\n';
  std::cout << "  -h,        --help               Print this help message."
            << '\n';
  std::cout << "  -p [path], --param-file [path]  Load parameter file." << '\n';
  std::cout << "  -c [path], --chkpt-file [path]  Load state from checkpoint."
            << '\n';
}

program_args::program_args(int argc, char **argv) {
  state parser_state = state::PROG_NAME;
  for (int idx = 0; idx < argc; ++idx) {
    auto arg = std::string(argv[idx]);
    switch (parser_state) {
    case state::PROG_NAME:
      prog = arg;
      parser_state = state::OPTION;
      break;
    case state::OPTION:
      if (arg == "-h" || arg == "--help")
        parser_state = state::HELP;
      else if (arg == "-p" || arg == "--param-file")
        parser_state = state::PARAM_PATH_VALUE;
      else if (arg == "-c" || arg == "--chkpt-file")
        parser_state = state::CHKPT_PATH_VALUE;
      else {
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
    case state::CHKPT_PATH_VALUE:
      if (!chkpt_path.has_value()) {
        chkpt_path = arg;
        parser_state = state::OPTION;
        break;
      } else {
        std::stringstream err_ss;
        err_ss << "Checkpoint path can be provided only once";
        error = err_ss.str();
        parser_state = state::INVALID;
      }
      break;
    case state::HELP:
    case state::INVALID:
      break;
    }
  }
}

void program::main(int argc, char **argv) {
  parse_args(argc, argv);
  setup_omp();

#pragma omp parallel default(none)
  {
    auto thr = thread(st);
    thr.main();
  }
}

void program::parse_args(int argc, char **argv) {
  auto args = program_args(argc, argv);
  if (args.error.has_value() || args.help) {
    args.print_help();
    auto exit_code = args.help ? EXIT_SUCCESS : EXIT_FAILURE;
    exit(exit_code);
  } else if (args.chkpt_path.has_value()) {
    throw std::runtime_error("Checkpoint are not (yet) supported.");
    exit(EXIT_FAILURE);
  } else {
    using namespace ioxx::xyaml;
    auto params_yml = defaults_yml();
    for (auto const &path : args.param_paths) {
      auto slice_yml = node::import(path);
      params_yml.merge(slice_yml);
    }
    params_yml >> st.params;
  }
}

void program::setup_omp() {
  omp_set_num_threads((int)st.params.gen.num_of_threads);
}
} // namespace cg::simul