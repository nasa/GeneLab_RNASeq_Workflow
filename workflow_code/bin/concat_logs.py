#! /usr/bin/env python
# Concatenate V&V logs

import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Concatenate V&V logs")
    parser.add_argument("--assay_suffix", default="", 
                        help="Assay suffix for output filename")
    args = parser.parse_args()

    logs = sorted(list(Path.cwd().glob("VV_in.csv*")))

    OUTPUT_FN = Path(f"VV_log_final{args.assay_suffix}.csv")

    for i, log in enumerate(logs):
        with open(log, "r") as in_f:
            if i == 0:  # first file
                contents = in_f.read()
            else:
                contents += "".join(in_f.readlines()[1:])  # skip header

    with open(OUTPUT_FN, "w") as out_f:
        out_f.write(contents)

if __name__ == "__main__":
    main()