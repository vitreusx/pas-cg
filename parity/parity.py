import argparse
import numpy as np
import pandas as pd
import subprocess
from dataclasses import dataclass
from pathlib import Path
from ruamel.yaml import safe_load


def parse_raw_data(raw_data: str):
    data = raw_data.split()
    idx = 0
    records = []

    while idx < len(data):
        traj_idx, step_idx, time, V = data[idx: idx + 4]
        traj_idx, step_idx, time, V = (
            int(traj_idx),
            int(step_idx),
            float(time),
            float(V),
        )
        records.append((traj_idx, step_idx, time, V))
        idx += 4

    df = pd.DataFrame.from_records(
        records, columns=["traj_idx", "step_idx", "time", "V"]
    )
    df.set_index(["traj_idx", "step_idx"])
    return df


def compare(tests):
    names = list(tests.keys())
    comp_df = pd.merge(
        tests[names[0]], tests[names[1]], how="outer",
        on=["traj_idx", "step_idx"]
    )
    comp_df[f"V_{names[0]}"] = comp_df["V_x"]
    comp_df[f"V_{names[1]}"] = comp_df["V_y"]
    comp_df["time"] = comp_df["time_x"]
    comp_df = comp_df.drop(["V_x", "V_y", "time_x", "time_y"], axis="columns")
    for name in names[2:]:
        comp_df = pd.merge(
            comp_df, tests[name][["traj_idx", "step_idx", "V"]], how="outer",
            on=["traj_idx", "step_idx"]
        )
        comp_df[f"V_{name}"] = comp_df["V"]
        comp_df = comp_df.drop(["V"], axis="columns")

    for name in names[1:]:
        v0 = comp_df[f"V_{names[0]}"]
        v1 = comp_df[f"V_{name}"]
        rel_diff = (v0 - v1).abs() / np.maximum(np.abs(v0), np.abs(v1))
        comp_df[f"rel_diff({names[0]},{name})"] = rel_diff

    return comp_df


@dataclass
class FortranOptions:
    num_of_threads: int
    variant: str

    @staticmethod
    def from_dict(data):
        return FortranOptions(
            num_of_threads=data["num of threads"],
            variant=data["variant"],
        )


def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--bin",
        required=True,
        help="path to the binary dir",
    )
    p.add_argument(
        "--test",
        required=True,
        help="path to test directory",
    )
    p.add_argument(
        "--source",
        required=True,
        help="path to original test directory, to put comp.csv to",
    )

    args = p.parse_args()
    bin = Path(args.bin).absolute()
    test = Path(args.test).absolute()
    src = Path(args.source).absolute()

    def run_cxx(input_file: Path):
        cxx = bin / "cg-cxx"
        subprocess.call([str(cxx), "-p", input_file.name],
                        cwd=input_file.parent)

    def run_fort(input_file: Path, suffix: str):
        options_path = input_file.parent / f"options{suffix}.yml"
        if not options_path.exists():
            options_path = Path(__file__).parent / "options.yml"

        with open(options_path, mode="r") as options_file:
            options = FortranOptions.from_dict(
                data=safe_load(options_file),
            )

        fort = bin / f"cg-fort-{options.variant}"

        if options.num_of_threads > 1:
            cmd = f"ulimit -s unlimited && env OMP_NUM_THREADS={options.num_of_threads} {str(fort)} {str(input_file.name)}"
            subprocess.call(cmd, shell=True, cwd=input_file.parent)
        else:
            cmd = [str(fort), str(input_file.name)]
            subprocess.call(cmd, cwd=input_file.parent)

    tests = {}
    for input_file in test.glob("**/inputfile*"):
        parent_dir_name = input_file.parent.name

        if parent_dir_name.startswith("cxx"):
            suffix = input_file.with_suffix("").name[len("inputfile"):]
            run_cxx(input_file)
        elif parent_dir_name.startswith("fort"):
            suffix = input_file.name[len("inputfile"):]
            run_fort(input_file, suffix)

        data_txt = open(input_file.parent / "raw_data.txt").read()
        test_data = parse_raw_data(data_txt)

        test_name = f"{parent_dir_name}{suffix}"
        tests[test_name] = test_data

    comp_df = compare(tests)
    comp_df.to_csv(test / "comp.csv", index=False, header=True)
    comp_df.to_csv(src / "comp.csv", index=False, header=True)


if __name__ == "__main__":
    main()
