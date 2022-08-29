# User Guide

## Introduction

## Installation

### Requirements

- Linux operating system;
- Git for cloning the repository, unless one uses a fixed release of the program;
- `gcc++` or `clang++` with ISO C++17 support;
- CMake 3.14+;
- (Optionally) Ninja build system;
- OpenMP;
- `xxd` program installed;
- Python 3, if using the [build script](build_script);
- Docker, if using the [Docker route](docker).

### Getting the program

To get the source code, you can clone the [project repository](https://github.com/vitreusx/pas-cg). For an end-user (i.e. just using the program and not extending it), shallow clone should be sufficient.

```shell
% git clone \
    --depth 1 \
    https://github.com/vitreusx/pas-cg.git

% cd pas-cg
```

### Compiling

(build_script)=
A build script has been provided, with following options:

```shell
% python3 build.py --help
usage: build.py [-h] [--target {cg,docs,parity,all,clean}] [--build-type {debug,release,staging}] [--generator {make,ninja}] [--single-file] [--use-mixed-precision] dir

Build script for PAS CG program

positional arguments:
  dir                   output directory

options:
  -h, --help            show this help message and exit
  --target {cg,docs,parity,all,clean}
                        CMake target to build (default: cg)
  --build-type {debug,release,staging}
                        build type (default: release)
  --generator {make,ninja}
                        build generator (default: ninja)
  --single-file         compile the program as a single object file, instead of separately compiling each source file and linking them (default: False)
  --use-mixed-precision
                        use floats for computing the forces and doubles for integration,instead of using doubles throughout (default: False)
```

(docker)=
### Docker

If all else fails, one can build and run the program as a Docker container. To build the image, run

```shell
$ ./docker/build.sh
```

in the root repo directory. To run the program, you can use

```shell
$ ./docker/run.sh [args]
```

with the arguments as usual. The run script prints the auto-generated container name - this can be used, for example, to kill the container via `docker rm -f ${name}`, as the Docker-run program cannot be stopped with the usual `Ctrl-c`.

## Quick Start

The program has (only) a command-line interface:

```shell
$ ./cg -h
PAS CG Coarse Grained Model
Documentation: https://vitreusx.github.io/pas-cg
Usage: ./cg [options...]
  -h,        --help               Print this help message.
  -p [path], --param-file [path]  Load parameter file.
  -c [path], --ckpt-file [path]   Load simulation from checkpoint.
```

One can either run the program from a checkpoint, or start with a *parameter file* (henceforth also called an *input file*). This file, in YAML format, completely specifies the simulation that will be run.

### Examples

To quickly check how the program functions, one can use provided example files, placed in the `cg/data` directory. For example, one can run

```shell
$ ./cg -p data/1ubq/inputfile.yml
```
to run a simulation of ubiquitin.

## Functionality

The simulation program can be described to function in a following fashion:

1. Arguments are read, and either the simulation is started from a [checkpoint](_ckpt), or an [input file](_inputfile) is read and run from that.
2. [The protein(s) are loaded](_input), either from a PDB file or a "sequence file" in a custom format.
3. The simulation proper is executed. The entire execution consists of a number of trajectories, with possibly different settings. Each trajectory starts from a (native or generated) conformation of the proteins loaded, which evolves according to the various force fields present and external factors for a given maximum simulation time.
   1. The peptide chains can be contained in a [simulation box](_sbox), bounded by a solid wall of potential. The geometry of the box can also be periodic.
   2. There are a number of force fields available:
      1. "Local" forces, i.e. ones which apply to residues separated by at most two other residues:
         1. [Tether forces](_tether) between connected bonds.
         2. [Chirality forces](_chir), if starting from a native conformation.
         3. [Bond angle forces](_angle), acting on triples of consecutive residues. The value can be either derived from the native structure, or a heurestic can be used if no such native value is available.
         4. [Dihedral angle forces](_angle), acting on quadruples of consecutive residues. Similarly as with the bond angle forces, there are native and heurestic versions of this potential.
         5. [Local repulsive forces](_lrep), acting on residues separated by a single residue.
      2. "Nonlocal" forces:
         1. [Pauli exclusion forces](_pauli), to prevent residues from passing through each other.
         2. [Debye-Hueckel electrostatic forces](_dh), with a constant-electric-permittivity and a relative-electric-permittivity versions.
         3. [Go model contacts](_nc), in which contacts in the native conformation act through the Lennard-Jones potential.
      3. [*Quasi-adiabatic potential*](_qa), in which native-like contacts are continuously formed and broken, according to a number of criteria.
      4. [*Pseudo-improper dihedral potential*](_pid), in which, depending on the geometry of the protein chains, the residues act as though there was a contact between them.
   3. There are a number of simulation types:
      1. "Free-form" simulations, in which the system undergoes evolution without any external intervention (except for the thermostat);
      2. [Simulations with Atomic Force Microscope (AFM)](_afm).
         1. (Optional) At the very start of the simulation, the terminal residue of the first chain present is pulled in the direction opposite to the direction linking it with the initial residue. This takes place for a specified amount of time, at which poiint the chain is released.
         2. The system undergoes *equilibration*, during which the only external forces acting upon the system are the solid walls of the simulation box, if present, and the thermostat.
         3. (Optional) After the equilibration, the initial and terminal residues of the first chain are pulled apart from each other. This continues until the forces become too large and the simulation ends prematurely.
      3. Simulations with deforming simulation box:
         1. (Optional) After the initial equilibration, the simulation box [is squeezed](_sbox#squeezing) to reach desired size or residue density. This is succeeded by another period of equilibration.
         2. (Optional) The simulation box is [relaxed until the forces acting on it are approximately zero](_sbox#finding-force-minimum). This is succeeded by another period of equilibration.
         3. (Optional) The system undergoes a series of [*oscillations*](_sbox#oscillations), during which the $z$-axis walls move towards each other, or sideways along the $x$ axis. This occurs with a given displacement amplitude and oscillation angular frequency for a number of cycles. Each such cycle is succeeded by a period of equilibration.
         4. (Optional) Then, at the end, [the simulation box is dragged in a single direction](_sbox#pulling-at-the-end).
         4. Throughout all this, at certain points the walls can become attractive, at which point the residues close to them become attracted and/or attached to them.
      4. [(Un)folding study](_nc#foldingunfolding-study):
         1. Unfolding study: We start from a native conformation and wait until all the native contacts have been broken.
         2. Folding study: We start from a random conformation and wait until all the native contacts have been restored.
         3. This can be repeated for a number of trajectories, and median (un)folding times can be measured.
4. [The output of the program](_output) consists of:
   1. A number of files about the current state of the system, emitted at regular intervals:
      1. An `*.out` file, which contains scalars such as the current time and energy;
      2. A `*.map` file with details concerning the contacts present and the values of bond/dihedral angles;
      3. A `*.pdb` file with current positions of the residues.
   2. (Optional) [A progress bar](_pbar).

To get started and run your own simulation, either:
- start from a given [checkpoint](_ckpt);
- write your own [input file](_inputfile).