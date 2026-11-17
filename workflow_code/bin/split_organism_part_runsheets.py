#!/usr/bin/env python3
"""Split a runsheet by Factor Value[organism part] for per-part DGE.

When that factor has 2+ values, write one runsheet per part and drop the
column so DGE treats it as a dataset split, not a grouping factor.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

ORGANISM_PART_NAMES = {"organism part"}


def organism_part_label(value: str) -> str:
    """Factor value for filenames: spaces become underscores."""
    return str(value).strip().replace(" ", "_")


def factor_value_columns(columns) -> list[str]:
    return [c for c in columns if str(c).startswith("Factor Value[") and str(c).endswith("]")]


def factor_name(col: str) -> str:
    return col[len("Factor Value[") : -1]


def find_organism_part_column(columns) -> str | None:
    for col in factor_value_columns(columns):
        if factor_name(col).strip().lower() in ORGANISM_PART_NAMES:
            return col
    return None


def unique_part_values(rows: list[dict], col: str) -> list[str]:
    vals = []
    seen = set()
    for row in rows:
        v = (row.get(col) or "").strip()
        if v not in seen:
            seen.add(v)
            vals.append(v)
    return vals


def part_label(value: str) -> str:
    return "_" + organism_part_label(value)


def _one_job(rows):
    return [{"output_label": "", "rows": rows, "drop_cols": []}]


def split_jobs(rows: list[dict], columns: list[str]) -> tuple[str | None, list[dict]]:
    col = find_organism_part_column(columns)
    if not col:
        return None, _one_job(rows)
    values = [v for v in unique_part_values(rows, col) if v]
    if len(values) < 2:
        return None, _one_job(rows)
    if any(not (row.get(col) or "").strip() for row in rows):
        raise ValueError(f"{col} contains blank values; assign an organism part to every sample before splitting DGE.")
    labels = [organism_part_label(v) for v in values]
    if len(set(labels)) != len(labels):
        raise ValueError(f"{col} values produce duplicate output labels after replacing spaces with underscores.")
    for value, label in zip(values, labels):
        if not re.fullmatch(r"[\w.-]+", label) or label == "rRNArm" or "_rRNArm" in label:
            raise ValueError(f"Unsupported organism-part output label for {value!r}; use letters, numbers, spaces, underscores, dots or hyphens, without the reserved 'rRNArm' label or '_rRNArm' suffix.")
    jobs = []
    for v in values:
        jobs.append({
            "output_label": part_label(v),
            "rows": [r for r in rows if (r.get(col) or "").strip() == v],
            "drop_cols": [col],
        })
    return factor_name(col), jobs


def write_runsheet(path: Path, rows: list[dict], columns: list[str], drop_cols: list[str]) -> None:
    keep = [c for c in columns if c not in drop_cols]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keep, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Split runsheet by organism part for DGE.")
    parser.add_argument("--runsheet", required=True, help="Input runsheet CSV")
    args = parser.parse_args()

    src = Path(args.runsheet)
    with src.open(newline="") as f:
        reader = csv.DictReader(f)
        columns = list(reader.fieldnames or [])
        rows = list(reader)

    factor, jobs = split_jobs(rows, columns)

    with open("stratify_by.txt", "w") as f:
        f.write(factor or "")

    with open("dge_jobs.tsv", "w") as manifest:
        manifest.write("output_label\tfilename\n")
        for job in jobs:
            label = job["output_label"]
            fname = f"dge_part{label or '_all'}.csv"
            write_runsheet(Path(fname), job["rows"], columns, job["drop_cols"])
            manifest.write(f"{label}\t{fname}\n")


if __name__ == "__main__":
    main()
