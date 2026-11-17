"""Shared raw/trimmed FASTQ gzip+format VV from per-file reports."""
from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path


def expected_fastqs(outdir, samples, paired_end, assay_suffix, kind):
    if kind == "raw":
        sub = ("00-RawData", "Fastq")
        suf = "raw"
    else:
        sub = ("01-TG_Preproc", "Fastq")
        suf = "trimmed"
    d = os.path.join(outdir, *sub)
    pe = paired_end[0] if len(paired_end) else False
    paths = []
    for sample in samples:
        if pe:
            paths.append(os.path.join(d, f"{sample}{assay_suffix}_R1_{suf}.fastq.gz"))
            paths.append(os.path.join(d, f"{sample}{assay_suffix}_R2_{suf}.fastq.gz"))
        else:
            paths.append(os.path.join(d, f"{sample}{assay_suffix}_{suf}.fastq.gz"))
    return paths


def load_reports(report_dir):
    p = Path(report_dir)
    files = [p] if p.is_file() else [f for f in p.rglob("*") if f.is_file()]
    reports = {}
    for f in files:
        for line in f.read_text(errors="replace").splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            reports[parts[0]] = {
                "gzip": parts[1],
                "format": parts[2],
                "msg": parts[3] if len(parts) > 3 else "",
            }
    return reports


def _script():
    p = Path(__file__).resolve().parent / "vv_check_fastq.sh"
    return str(p) if p.exists() else "vv_check_fastq.sh"


def collect_reports(paths, report_dir=None):
    if report_dir:
        return load_reports(report_dir)
    tmp = tempfile.mkdtemp(prefix="vv_fq_")
    script = _script()
    for path in paths:
        if not os.path.exists(path):
            continue
        out = Path(tmp) / f"{os.path.basename(path)}.vv_fastq.tsv"
        with open(out, "w") as fh:
            subprocess.run([script, path], stdout=fh, check=False)
    return load_reports(tmp)


def apply_fastq_check_reports(log_path, component, expected_paths, reports, log_check_result):
    checked = [p for p in expected_paths if os.path.exists(p)]
    if not checked:
        log_check_result(
            log_path, component, "all", "check_gzip_integrity", "HALT",
            "No FASTQ files found to check", "",
        )
        log_check_result(
            log_path, component, "all", "validate_fastq_format", "HALT",
            "No files found to validate" if component == "raw_reads"
            else "No trimmed FASTQ files found to validate",
            "",
        )
        return False

    gzip_fail = []
    fmt_fail = []
    missing = []
    for path in checked:
        r = reports.get(os.path.basename(path))
        if r is None:
            missing.append(path)
            continue
        if r["gzip"] != "OK":
            gzip_fail.append(path)
        if r["format"] != "OK":
            fmt_fail.append(path)

    gzip_fail.extend(missing)
    fmt_fail.extend(missing)
    n = len(checked)
    ok = True
    if gzip_fail:
        print(f"GZIP integrity check failed for {len(gzip_fail)} of {n} files")
        for f in gzip_fail:
            print(f"  - {f}")
        log_check_result(
            log_path, component, "all", "check_gzip_integrity", "HALT",
            f"Integrity check failed for {len(gzip_fail)} of {n} files",
            ",".join(gzip_fail),
        )
        ok = False
    else:
        print(f"GZIP integrity check passed for all {n} FASTQ files")
        log_check_result(
            log_path, component, "all", "check_gzip_integrity", "GREEN",
            f"Integrity check passed for all {n} files", "",
        )

    if fmt_fail:
        print(f"WARNING: The following FASTQ files are invalid:")
        for f in fmt_fail:
            print(f"  - {f}")
        log_check_result(
            log_path, component, "all", "validate_fastq_format", "HALT",
            f"Invalid format in {len(fmt_fail)} of {n} files",
            ",".join(fmt_fail),
        )
        ok = False
    else:
        print(f"All {n} FASTQ files are valid")
        log_check_result(
            log_path, component, "all", "validate_fastq_format", "GREEN",
            "All files valid", "",
        )
    return ok


def check_fastq_gzip_and_format(
    outdir, samples, paired_end, log_path, assay_suffix, component, kind,
    report_dir, log_check_result,
):
    paths = expected_fastqs(outdir, samples, paired_end, assay_suffix, kind)
    reports = collect_reports(paths, report_dir)
    return apply_fastq_check_reports(
        log_path, component, paths, reports, log_check_result,
    )
