import argparse
import os
import shutil
import subprocess as sp
import tempfile
from pathlib import Path


def make_parser():
    p = argparse.ArgumentParser(
        prog="build.py",
        description="Build script for PAS CG program",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    mode = p.add_subparsers(dest="mode")

    cxx = mode.add_parser("cxx")

    cxx.add_argument(
        "-t", "--target",
        choices=["cg", "docs", "parity", "all"],
        default="cg",
        help="CMake target to build",
    )
    cxx.add_argument(
        "-B", "--build-type",
        choices=["debug", "release", "staging"],
        default="release",
        help="build type",
    )
    cxx.add_argument(
        "-G", "--generator",
        choices=["make", "ninja"],
        default="ninja",
        help="build generator",
    )
    cxx.add_argument(
        "-C", "--compiler",
        choices=["clang", "gcc", "intel"],
        default="clang",
        help="compiler choice"
    )
    cxx.add_argument(
        "--single-file",
        default=False,
        help="whether to compile the program as a single object file, instead of "
             "separately compiling each source file and linking them",
    )
    cxx.add_argument(
        "--mixed-precision",
        default=False,
        help="whether to  use floats for computing the forces and doubles for integration,"
             "instead of using doubles throughout",
    )
    cxx.add_argument(
        "--vectorized-impls",
        default=True,
        help="whether to use vectorized implementations of the algorithms",
    )

    cxx.add_argument("-f", "--output-file",
                     help="output executable file (cg)")

    cxx.add_argument("-d", "--output-dir",
                     help="output CMake directory")

    fort = mode.add_parser("fort")

    fort.add_argument(
        "-B", "--build-type",
        choices=["debug", "release-g", "release"],
        default="release",
        help="build type",
    )

    fort.add_argument(
        "-C", "--compiler",
        choices=["gcc", "intel"],
        default="gcc",
        help="compiler to use",
    )

    fort.add_argument(
        "--variant",
        choices=["pho5", "pho6"],
        default="pho6",
        help="Fortran program variant",
    )

    fort.add_argument(
        "output_file",
        default="cg",
        help="executable path",
    )

    return p


def format_cmd(cmd) -> str:
    def enclose(s):
        return f'"{s}"' if " " in s else s

    cmd = [cmd[0], *map(enclose, cmd[1:])]
    cmd = " ".join(cmd)
    cmd = f"% {cmd}"
    return cmd


def execute(cmd):
    print(format_cmd(cmd))

    code = sp.call(cmd)
    if code != 0:
        print(f"Command {cmd} failed with code {code}")
        exit(1)


def to_bool(s):
    return int(s in ("1", "true", "True", "ON", True))


def build_cxx(args):
    root_dir = Path(__file__).parent

    if args.output_dir is not None:
        out_dir = Path(args.output_dir)
        print(format_cmd(["mkdir", "-p", str(out_dir)]))
        out_dir.mkdir(exist_ok=True, parents=True)
    else:
        tmp_out_dir = tempfile.TemporaryDirectory()
        out_dir = Path(tmp_out_dir.name)

    if args.output_file is not None:
        out_file = Path(args.output_file)

    cmd = ["cmake"]

    if args.generator == "ninja" and shutil.which("ninja") is not None:
        generator = "Ninja"
    else:
        generator = "Unix Makefiles"
    cmd += ["-G", generator]

    if args.compiler == "clang" and shutil.which("clang") is not None:
        cmd += ["-DCMAKE_C_COMPILER=clang", "-DCMAKE_CXX_COMPILER=clang++"]
    elif args.compiler == "intel" and shutil.which("icc") is not None:
        cmd += ["-DCMAKE_C_COMPILER=icc", "-DCMAKE_CXX_COMPILER=icc"]
    elif shutil.which("gcc") is not None:
        cmd += ["-DCMAKE_C_COMPILER=gcc", "-DCMAKE_CXX_COMPILER=g++"]
    else:
        raise RuntimeError("Could not find any compiler!")

    cmd += ["-B", str(out_dir)]
    cmd += ["-S", str(root_dir)]

    build_type_map = {"debug": "Debug", "release": "Release",
                      "staging": "Staging"}
    build_type = build_type_map[args.build_type]
    cmd += [f"-DCMAKE_BUILD_TYPE={build_type}"]

    cmd += [
        f"-DSINGLE_FILE={to_bool(args.single_file)}",
        f"-DUSE_MIXED_PRECISION={to_bool(args.mixed_precision)}",
        f"-DUSE_VECTORIZED_IMPLS={to_bool(args.vectorized_impls)}"
    ]

    execute(cmd)

    cmd = ["cmake"]
    cmd += ["--build", str(out_dir)]
    cmd += ["--parallel", str(os.cpu_count())]
    cmd += ["--target", args.target]

    execute(cmd)

    if args.output_file is not None:
        shutil.copy2(out_dir / "cg" / "cg", out_file)

    if args.output_dir is None:
        tmp_out_dir.cleanup()


def build_fort(args):
    cmd = []

    if args.compiler == "intel" and shutil.which("ifort") is not None:
        cmd += ["ifort", "-std=legacy", "-mcmodel=large", "-qopenmp"]
        if args.build_type != "debug":
            cmd += ["-fast"]
    elif args.compiler == "gcc" and shutil.which("gfortran") is not None:
        cmd += ["gfortran", "-std=legacy", "-mcmodel=large", "-fopenmp"]
        if args.build_type != "debug":
            cmd += ["-Ofast", "-march=native"]
    else:
        raise RuntimeError("Could not find compiler.")

    if args.build_type != "release":
        cmd += ["-g"]

    out_f = Path(args.output_file)
    cmd += ["-o", str(out_f)]

    root_dir = Path(__file__).parent
    cg_f = root_dir / "fort" / f"cg_{args.variant}.f"
    cmd += [str(cg_f)]

    execute(cmd)


def main():
    args = make_parser().parse_args()
    if args.mode == "cxx":
        build_cxx(args)
    elif args.mode == "fort":
        build_fort(args)


if __name__ == "__main__":
    main()
