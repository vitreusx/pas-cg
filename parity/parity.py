import argparse
import pandas as pd
import subprocess
from pathlib import Path


def parse_raw_data(raw_data: str):
    data = raw_data.split()
    idx = 0
    records = []

    while idx < len(data):
        traj_idx, step_idx, V = data[idx: idx + 3]
        traj_idx, step_idx, V = int(traj_idx), int(step_idx), float(V)
        records.append((traj_idx, step_idx, V))
        idx += 3

    df = pd.DataFrame.from_records(records,
                                   columns=["traj_idx", "step_idx", "V"])
    df.set_index(["traj_idx", "step_idx"])
    return df


def compare(test, fort_data, cxx_data):
    comp_df = pd.merge(cxx_data, fort_data, how="outer",
                       on=["traj_idx", "step_idx"])
    comp_df["V_cxx"] = comp_df["V_x"]
    comp_df["V_fort"] = comp_df["V_y"]
    comp_df = comp_df.drop(["V_x", "V_y"], axis="columns")
    comp_df["rel_diff"] = (comp_df["V_cxx"] - comp_df["V_fort"]).abs() / \
                          comp_df["V_cxx"].abs()
    return comp_df


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-c", "--cxx", required=True, help="path to cxx program")
    p.add_argument("-f", "--fort", required=True,
                   help="path to fortran program")
    p.add_argument("-t", "--test", required=True, help="path to test directory")
    p.add_argument("-s", "--source", required=True,
                   help="path to original test directory, to put comp.csv to")

    args = p.parse_args()
    cxx = Path(args.cxx)
    fort = Path(args.fort)
    test = Path(args.test)
    src = Path(args.source)

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

    comp_df = compare(test, fort_data, cxx_data)
    comp_df.to_csv(test / "comp.csv", index=False, header=True)
    comp_df.to_csv(src / "comp.csv", index=False, header=True)
    

if __name__ == "__main__":
    main()
