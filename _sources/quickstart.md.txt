# Quick Start

## Requirements

- Linux operating system;
- Git;
- `gcc++` or `clang++` with C++17 support;
- CMake 3.14+;
- OpenMP;
- `xxd` program installed.

## Getting the program

### Releases

For convenience, a number of releases are available at
the [Github page](https://github.com/vitreusx/pas-cg/releases). The releases
contain:

- the executable `cg`, compiled on the latest Ubuntu version;
- `data` directory with examples.

For more info, see
the [workflow file](https://github.com/vitreusx/pas-cg/blob/main/.github/workflows/create-release.yml)
.

```{warning}
Some shared libraries may not have the same name as on the Ubuntu machine, resulting in an error. To fix it, check the offending libraries with `ldd cg`.
```

### Compiling from source

#### Getting the source code

To get the source code, clone
the [vitreusx/pas-cg](https://github.com/vitreusx/pas-cg) repository and pull
the submodules. For regular use, shallow clones are sufficient:

```shell
git clone --depth 1 \
  --recurse-submodules \
  --shallow-submodules \
  https://github.com/vitreusx/pas-cg.git
cd pas-cg/
```

#### Compiling

The project is built using [CMake](https://cmake.org/). The two targets of
interest are:

- `cg`, the program proper;
- `docs`, target for generating the documentation.

The `cg` target supports following configurations:

- `Debug`: version without optimizations and with debug info attached;
- `Release`: fully optimized version;
- `Staging`: `Release` with debug info attached;
- `Deployment`: version with optimizations but without `-march=native`.

For example, to compile the program in the Release mode with `make` to
the `release/` directory, enter:

```shell
mkdir release
cd release
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4 cg
```

## Running

```
$ ./cg --help
PAS CG Coarse Grained Model
Documentation: https://vitreusx.github.io/pas-cg
Usage: /home/talos/next/code/pas-cg/staging/cg/cg [options...]
  -h,        --help               Print this help message.
  -p [path], --param-file [path]  Load parameter file.
  -c [path], --chkpt-file [path]  Load state from checkpoint.
  -o [path], --output     [path]  Set the output directory.
```

The input to the program consist of the parameter file(s), which allow the user
to specify the various parameters of the simulation. The (optional) output
consists of a progress bar and files in the output directory. For example, the
command to run an example simulation could look like:

```shell
./cg -p data/${example}/inputfile.yml
```

The provided parameter/input files contain only the differences from the input
file with the defaults, which is included by default for the sake of convenience
and is located
at [`cg/defaults.yml`](https://github.com/vitreusx/pas-cg/blob/main/cg/default.yml)
. To be precise, the implicit default parameter file gets overriden by the
user-provided parameter files in the order in which they are listed. The output
directory can be overriden with the `-o` option.

## Examples

A number of examples are provided for the purposes of demonstration in
the [`cg/data/`](https://github.com/vitreusx/pas-cg/tree/main/cg/data)
directory.

- [`1ubq`](https://github.com/vitreusx/pas-cg/tree/main/cg/data/1ubq): folding
  of ubiquitin, computation of median required time for re-establishing the
  native contacts;
- [`9aac`](https://github.com/vitreusx/pas-cg/tree/main/cg/data/9aac):
  simulation of Atomic Force Microscope (AFM) stretching of an intrinsically
  disordered protein, $\alpha$-synuclein, with the pseudo-improper dihedral
  custom potential;
- [`glut`](https://github.com/vitreusx/pas-cg/tree/main/cg/data/1ubq):
  simulation of deforming a box filled with partially structured gluten
  proteins, with the quasi-adiabatic custom potential.

```{warning}
`glut` example is, at the moment, unfinished.
```

## Docker

```{warning}
This method is somewhat experimental.
```

"As a last resort", one can build and run the program as a Docker container. To
build the image, run the following commands in the `pas-cg/` directory:

```shell
docker build --tag "vitreusx/pas-cg:latest" .
```

To run the program, execute:

```shell
docker run -it \
  --name ${name} \
  -v $(pwd):/host \
  --user "$(id -u):$(id -g)" \
  "vitreusx/pas-cg:latest" \
  [args...]
```

where `${name}` is (optional) human-readable name for the container. These two
commands are provided as `docker-build` and `docker-run` scripts in the root
directory. To run them, Bash is required. `docker-build` must be run in the root
directory. The usage of `docker-run` is:

```
Usage:
./docker-run name [args...]
```

```{warning}
The container cannot be stopped with Ctrl-c. To kill it, use `docker rm -f ${name}`.
```