# Quick Start

## Requirements

- Linux operating system;
- `gcc++` or `clang++` with C++17 support;
- CMake 3.14+;
- OpenMP.

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
consists of a progress bar and files in the output directory. For example, the
command to run an example simulation could look like:

```shell
./cg data/${example}/inputfile.yml
```

The provided parameter/input files contain only the differences from the input
file with the defaults, which is included by default for the sake of convenience
and is located at `cg/src/simul/defaults.yml`. To be precise, the implicit
default parameter file gets overriden by the user-provided parameter files in
the order in which they are listed.

## Examples

A number of examples are provided for the purposes of demonstration:

- `1ubq`: folding of ubiquitin, computation of median required time for
  re-establishing the native contacts;
- `9aac`: simulation of Atomic Force Microscope (AFM) stretching of an
  intrinsically disordered protein, $\alpha$-synuclein, with the pseudo-improper
  dihedral custom potential;
- `glut`: simulation of deforming a box filled with partially structured gluten
  proteins, with the quasi-adiabatic custom potential.

```{warning}
`glut` example is, at the moment, unfinished.
```