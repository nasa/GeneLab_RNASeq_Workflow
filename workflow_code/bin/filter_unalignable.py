#!/usr/bin/env python3
"""Drop runsheet rows whose RSEM pct_unalignable is at or above a threshold."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


def sample_name_column(columns: list[str]) -> str:
    for c in columns:
        if str(c).lower() in ("sample name", "sample.name"):
            return c
    for c in columns:
        cl = str(c).lower()
        if "sample" in cl and "name" in cl:
            return c
    raise SystemExit("Runsheet has no Sample Name column")


def pct_unalignable_from_mqc(root: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    if not root.exists():
        return out
    for path in root.rglob("multiqc_data.json"):
        try:
            raw = json.loads(path.read_text()).get("report_saved_raw_data") or {}
            rsem = raw.get("multiqc_rsem") or {}
        except (OSError, json.JSONDecodeError, AttributeError):
            continue
        for sample, counts in rsem.items():
            try:
                total = (
                    counts["Unique"]
                    + counts["Multi"]
                    + counts["Filtered"]
                    + counts["Unalignable"]
                )
                if total <= 0:
                    continue
            except (KeyError, TypeError):
                continue
            out[str(sample)] = counts["Unalignable"] / total * 100
    return out


def strip_suffix(name: str, assay_suffix: str) -> str:
    name = str(name).strip()
    if assay_suffix and name.endswith(assay_suffix):
        return name[: -len(assay_suffix)]
    return name


def exclude_names(runsheet_names: list[str], pct: dict[str, float], threshold: float, assay_suffix: str) -> list[tuple[str, float]]:
    by_stripped = {strip_suffix(k, assay_suffix): v for k, v in pct.items()}
    dropped = []
    for name in runsheet_names:
        key = strip_suffix(name, assay_suffix)
        if key in by_stripped and by_stripped[key] >= threshold:
            dropped.append((name, by_stripped[key]))
        elif name in pct and pct[name] >= threshold:
            dropped.append((name, pct[name]))
    return dropped


def write_runsheet(path: Path, rows: list[dict], columns: list[str]) -> None:
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--runsheet", required=True)
    p.add_argument("--multiqc-data", default="", help="RSEM MultiQC data directory")
    p.add_argument("--threshold", type=float, default=60)
    p.add_argument("--assay-suffix", default="")
    args = p.parse_args()

    src = Path(args.runsheet)
    with src.open(newline="") as f:
        reader = csv.DictReader(f)
        columns = list(reader.fieldnames or [])
        rows = list(reader)
    sample_col = sample_name_column(columns)
    names = [(r.get(sample_col) or "").strip() for r in rows]

    mqc = Path(args.multiqc_data) if args.multiqc_data else Path()
    pct = pct_unalignable_from_mqc(mqc) if args.multiqc_data else {}
    if args.multiqc_data and not pct:
        print("WARNING: no RSEM unalignable rates in MultiQC data; keeping all samples", file=sys.stderr)

    dropped = exclude_names(names, pct, args.threshold, args.assay_suffix)
    drop_set = {n for n, _ in dropped}
    kept = [r for r in rows if (r.get(sample_col) or "").strip() not in drop_set]
    if not kept:
        raise SystemExit("All samples were at or above the unalignable threshold")

    write_runsheet(Path("dge_runsheet.csv"), kept, columns)
    with open("dropped_unalignable.tsv", "w") as f:
        f.write("sample\tpct_unalignable\n")
        for name, val in dropped:
            f.write(f"{name}\t{val:.4f}\n")

    if dropped:
        print(f"Dropped {len(dropped)} sample(s) at pct_unalignable >= {args.threshold:g}: {[n for n, _ in dropped]}")
    else:
        print("No samples dropped")


if __name__ == "__main__":
    main()
