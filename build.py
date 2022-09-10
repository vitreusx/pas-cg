import argparse
import os
import shutil
import subprocess as sp
from pathlib import Path


def make_parser():
    parser = argparse.ArgumentParser(
        prog="build.py",
        description="Build script for PAS CG program",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--target",
        choices=["cg", "docs", "parity", "all", "clean"],
        default="cg",
        help="CMake target to build",
    )
    parser.add_argument(
        "--build-type",
        choices=["debug", "release", "staging"],
        default="release",
        help="build type",
    )
    parser.add_argument(
        "--generator",
        choices=["make", "ninja"],
        default="ninja",
        help="build generator",
    )
    parser.add_argument(
        "--single-file",
        type=bool,
        default=false,
        help="compile the program as a single object file, instead of "
             "separately compiling each source file and linking them",
    )
    parser.add_argument(
        "--use-mixed-precision",
        type=bool,
        default=false,
        help="use floats for computing the forces and doubles for integration,"
             "instead of using doubles throughout",
    )
    parser.add_argument(
        "--use-vectorized-impls",
        type=bool,
        default=true,
        help="use vectorized implementations of the algorithms"
    )
    parser.add_argument("dir", default="build", help="output directory")

    return parser


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


def main():
    args = make_parser().parse_args()

    root_dir = Path(__file__).parent
    out_dir = Path(args.dir)

    print(format_cmd(["mkdir", "-p", str(out_dir)]))
    out_dir.mkdir(exist_ok=True, parents=True)

    cmd = ["cmake"]

    if args.generator == "ninja" and shutil.which("ninja") is not None:
        generator = "Ninja"
    else:
        generator = "Unix Makefiles"
    cmd += ["-G", generator]

    cmd += ["-B", str(out_dir)]
    cmd += ["-S", str(root_dir)]

    build_type_map = {"debug": "Debug", "release": "Release",
                      "staging": "Staging"}
    build_type = build_type_map[args.build_type]
    cmd += [f"-DCMAKE_BUILD_TYPE={build_type}"]

    cmd += [
        f"-DSINGLE_FILE={int(args.single_file)}",
        f"-DUSE_MIXED_PRECISION={int(args.use_mixed_precision)}",
        f"-DUSE_VECTORIZED_IMPLS={int(args.use_vectorized_impls)}"
    ]

    execute(cmd)

    cmd = ["cmake"]
    cmd += ["--build", str(out_dir)]
    cmd += ["--parallel", str(os.cpu_count())]
    cmd += ["--target", args.target]

    execute(cmd)


if __name__ == "__main__":
    main()
