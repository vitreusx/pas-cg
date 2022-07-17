import argparse
from dataclasses import dataclass
from pathlib import Path
import subprocess
import tempfile
import shutil
import numpy as np


@dataclass
class Step:
    V: float
    r: np.ndarray
    F: np.ndarray


def parse_raw_data(raw_data: str):
    data = raw_data.split()
    idx = 0
    result = {}
    while idx < len(data):
        traj_idx, step_idx, V, num_res = data[idx : idx + 4]
        idx += 4
        step_idx, V, num_res = int(step_idx), float(V), int(num_res)
        r, F = np.zeros((num_res, 3)), np.zeros((num_res, 3))
        for res_idx in range(num_res):
            rx, ry, rz, Fx, Fy, Fz = (float(x) for x in data[idx : idx + 6])
            idx += 6
            r[res_idx] = np.array([rx, ry, rz])
            F[res_idx] = np.array([Fx, Fy, Fz])
        result[traj_idx, step_idx] = Step(V, r, F)
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-c", "--cxx", required=True, help="path to cxx program")
    p.add_argument("-f", "--fort", required=True, help="path to fortran program")
    p.add_argument("-t", "--test", required=True, help="path to test directory")

    args = p.parse_args()
    cxx = Path(args.cxx)
    fort = Path(args.fort)
    test = Path(args.test)

    # with tempfile.TemporaryDirectory() as tmpdir:
    #     tmpdir = Path(tmpdir)
    #     shutil.copytree(test, tmpdir, dirs_exist_ok=True)

    #     subprocess.call([cxx, "-p", "inputfile.yml"], cwd=tmpdir / "cxx")
    #     cxx_data_txt = open(tmpdir / "cxx" / "raw_data.txt").read()
    #     cxx_data = parse_raw_data(cxx_data_txt)

    #     subprocess.call([fort, "inputfile"], cwd=tmpdir / "fort")
    #     fort_data_txt = open(tmpdir / "fort" / "raw_data.txt").read()
    #     fort_data = parse_raw_data(fort_data_txt)

    subprocess.call([cxx, "-p", "inputfile.yml"], cwd=test / "cxx")
    cxx_data_txt = open(test / "cxx" / "raw_data.txt").read()
    cxx_data = parse_raw_data(cxx_data_txt)

    subprocess.call([fort, "inputfile"], cwd=test / "fort")
    fort_data_txt = open(test / "fort" / "raw_data.txt").read()
    fort_data = parse_raw_data(fort_data_txt)


if __name__ == "__main__":
    main()
