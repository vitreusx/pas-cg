# Quick Start

## Requirements

- Linux operating system;
- `gcc++` or `clang++` with C++17 support;
- CMake 3.14+;
- OpenMP;
- Boost (more specifically the `program_options` component).

## Installation

To get the source code, clone
the [vitreusx/pas-cg](https://github.com/vitreusx/pas-cg) repository and pull
the submodules. For regular use, shallow clones are sufficient:

```shell
git clone --depth 1 --recurse-submodules --shallow-submodules https://github.com/vitreusx/pas-cg.git
cd pas-cg/
```

## Compiling

The project is built using [CMake](https://cmake.org/). The two targets of
interest are `cg`, the program proper, and `docs` for generating the
documentation. The `cg` target supports `Debug` and `Release` configurations.

```shell
mkdir release
cd release
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4 cg
```

## Running

```
Usage: cg [options] param-files...
Allowed options:
  --help                print this help message
  --param-files arg     parameter files
```

The input to the program consist of the parameter file(s), which allow the user
to specify the various parameters of the simulation. The (optional) output
consists of a progress bar and files in the output directory. Usually, there
will be two parameter files: the defaults (provided with the code), and the
overrides, so for example we may call:

```shell
./cg data/default/inputfile.yml data/${example}/inputfile.yml
```

```{warning}
All the parameters must be included. In particular, unless one modifies the 
defaults file directly, it is necessary to include it as well as the overrides file.
```

## Examples

A number of examples are provided for the purposes of demonstration:

- `1ubq`: folding of ubiquitin, computation of median required time for
  re-establishing the native contacts;
- `9aac`: simulation of Atomic Force Microscope (AFM) stretching of an
  intrinsically disordered protein, $\alpha$-synuclein, with the pseudo-improper
  dihedral custom potential;
- `glut`: simulation of deforming a box filled with partially structured gluten
  proteins, with the quasi-adiabatic custom potential.