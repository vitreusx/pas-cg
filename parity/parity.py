import argparse
import pandas as pd
import subprocess
from pathlib import Path
import numpy as np


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


def compare(tests):
    names = list(tests.keys())
    comp_df = pd.merge(tests[names[0]], tests[names[1]],
                       how="outer", on=["traj_idx", "step_idx"])
    comp_df[f"V_{names[0]}"] = comp_df["V_x"]
    comp_df[f"V_{names[1]}"] = comp_df["V_y"]
    comp_df = comp_df.drop(["V_x", "V_y"], axis="columns")
    for name in names[2:]:
        comp_df = pd.merge(comp_df, tests[name],
                           how="outer", on=["traj_idx", "step_idx"])
        comp_df[f"V_{name}"] = comp_df["V"]
        comp_df = comp_df.drop(["V"], axis="columns")

    for name in names[1:]:
        v0 = comp_df[f"V_{names[0]}"]
        v1 = comp_df[f"V_{name}"]
        rel_diff = (v0 - v1).abs() / np.maximum(np.abs(v0), np.abs(v1))
        comp_df[f"rel_diff({names[0]},{name})"] = rel_diff
    
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

    tests = {}
    for input_file in test.glob("**/inputfile*"):        
        parent_dir_name = input_file.parent.name
        if parent_dir_name.startswith("cxx"):
            subprocess.call([cxx, "-p", input_file.name], cwd=input_file.parent)
        elif parent_dir_name.startswith("fort"):
            subprocess.call([fort, "inputfile"], cwd=input_file.parent)
        
        data_txt = open(input_file.parent / "raw_data.txt").read()
        test_data = parse_raw_data(data_txt)

        input_file_suffix = input_file.with_suffix("").name[len("inputfile"):]
        test_name = f"{parent_dir_name}{input_file_suffix}"
        tests[test_name] = test_data

    comp_df = compare(tests)
    comp_df.to_csv(test / "comp.csv", index=False, header=True)
    comp_df.to_csv(src / "comp.csv", index=False, header=True)
    

if __name__ == "__main__":
    main()
