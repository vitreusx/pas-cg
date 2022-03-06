# Tutorial

## Requirements

- Linux operating system
- `gcc++` or `clang++` with C++17 support
- CMake 3.14+
- OpenMP

## Installation

To get the source code, clone the `https://github.com/vitreusx/pas-cg` repository and pull the submodules. For regular use, shallow clones are sufficient:

```zsh
git clone --depth 1 --recurse-submodules --shallow-submodules https://github.com/vitreusx/pas-cg.git
cd pas-cg/
```

## Compiling

The project is built using CMake. The two targets of interest are `cg`, the program proper, and `docs` for generating the documentation. The `cg` target supports `Debug` and `Release` configurations.

```zsh
mkdir release
cd release
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4 cg
```
