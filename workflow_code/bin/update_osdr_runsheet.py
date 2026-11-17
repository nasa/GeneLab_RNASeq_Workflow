#!/usr/bin/env python3

import argparse
import csv
import sys
from pathlib import Path


def load_file_list(path):
    mapping = {}
    with path.open(newline="") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line.strip():
                continue
            name, tab, rest = line.partition("\t")
            if not tab or name == "Filename":
                continue
            url = rest.split("\t", 1)[0]
            if name in mapping and mapping[name] != url:
                _die(f"Multiple URLs found for filename {name}")
            mapping[name] = url
    return mapping


def _truthy(value):
    return str(value).strip().lower() == "true"


def _die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    raise SystemExit(1)


def suffixes(suffix):
    return [suffix, ""] if suffix else [""]


def resolve_trimmed(mapping, glds, sample, paired, suffix):
    for candidate in suffixes(suffix):
        prefix = f"{glds}_rna_seq_{sample}{candidate}"
        if paired:
            r1 = mapping.get(f"{prefix}_R1_trimmed.fastq.gz")
            r2 = mapping.get(f"{prefix}_R2_trimmed.fastq.gz")
            if r1 and r2:
                return {"trimmed_read1_path": r1, "trimmed_read2_path": r2}
        else:
            se = mapping.get(f"{prefix}_trimmed.fastq.gz") or mapping.get(f"{prefix}_R1_trimmed.fastq.gz")
            if se:
                return {"trimmed_read1_path": se, "trimmed_read2_path": ""}
    _die(f"Could not find trimmed reads for {sample} in file list")


def resolve_bam(mapping, glds, sample, mode, suffix):
    endings = ["_sorted.bam", ".bam"] if mode == "microbes" else ["_Aligned.toTranscriptome.out.bam"]
    for candidate in suffixes(suffix):
        for ending in endings:
            url = mapping.get(f"{glds}_rna_seq_{sample}{candidate}{ending}")
            if url:
                return {"bam_path": url}
    _die(f"Could not find BAM file for {sample} in file list")


def resolve_bam_pair(mapping, glds, sample, suffix):
    pairs = []
    for candidate in suffixes(suffix):
        prefix = f"{glds}_rna_seq_{sample}{candidate}"
        counting = mapping.get(prefix + "_Aligned.toTranscriptome.out.bam")
        for ending in ("_Aligned.sortedByCoord_sorted.out.bam", "_Aligned.sortedByCoord.out.bam", "_Aligned.out.bam"):
            genome = mapping.get(prefix + ending)
            if counting and genome:
                pairs.append({"bam_path": counting, "genome_bam_path": genome})
    if len(pairs) != 1:
        _die(f"BAM auto strandedness for {sample} needs one unambiguous transcriptome/genome pair; found {len(pairs)}. Supply paired runsheet paths or explicit --strandedness.")
    return pairs[0]


def resolve_genes(mapping, glds, sample, suffix):
    for candidate in suffixes(suffix):
        url = mapping.get(f"{glds}_rna_seq_{sample}{candidate}.genes.results")
        if url:
            return {"genes_results_path": url}
    _die(f"Could not find genes.results file for {sample} in file list")


def resolve_counts(mapping, glds, mode, suffix):
    counts_type = "FeatureCounts" if mode == "microbes" else "RSEM"
    for candidate in suffixes(suffix):
        url = mapping.get(f"{glds}_rna_seq_{counts_type}_Unnormalized_Counts{candidate}.csv")
        if url:
            return url
    _die(f"Could not find {counts_type} raw counts table in file list")


def resolve_dge(mapping, glds, suffix):
    for candidate in suffixes(suffix):
        url = mapping.get(f"{glds}_rna_seq_differential_expression{candidate}.csv")
        if url:
            return url
    _die("Could not find differential expression table in file list")


def update_runsheet(rows, fieldnames, *, mapping, glds, entry_point, mode, assay_suffix, strandedness="none"):
    extra = []
    if entry_point == "trimmed_reads":
        extra = ["trimmed_read1_path", "trimmed_read2_path"]
        for row in rows:
            sample = row["Sample Name"]
            row.update(resolve_trimmed(mapping, glds, sample, _truthy(row.get("paired_end")), assay_suffix))
    elif entry_point == "bam_files":
        extra = ["bam_path"]
        if mode != "microbes" and strandedness == "auto":
            extra.append("genome_bam_path")
        for row in rows:
            if mode != "microbes" and strandedness == "auto":
                row.update(resolve_bam_pair(mapping, glds, row["Sample Name"], assay_suffix))
            else:
                row.update(resolve_bam(mapping, glds, row["Sample Name"], mode, assay_suffix))
    elif entry_point == "genes_results":
        extra = ["genes_results_path"]
        for row in rows:
            row.update(resolve_genes(mapping, glds, row["Sample Name"], assay_suffix))
    elif entry_point == "counts_table":
        extra = ["counts_table_path"]
        url = resolve_counts(mapping, glds, mode, assay_suffix)
        for row in rows:
            row["counts_table_path"] = url
    elif entry_point == "dge_table":
        extra = ["dge_table_path"]
        url = resolve_dge(mapping, glds, assay_suffix)
        for row in rows:
            row["dge_table_path"] = url

    return rows, list(fieldnames) + [c for c in extra if c not in fieldnames]


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--runsheet", required=True, type=Path)
    p.add_argument("--file-list", required=True, type=Path)
    p.add_argument("--glds", required=True)
    p.add_argument("--entry-point", required=True, choices=["trimmed_reads", "bam_files", "genes_results", "counts_table", "dge_table"])
    p.add_argument("--mode", default="default", choices=["default", "microbes"])
    p.add_argument("--assay-suffix", default="")
    p.add_argument("--strandedness", default="none", choices=["auto", "none", "forward", "reverse"])
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args(argv)

    mapping = load_file_list(args.file_list)
    with args.runsheet.open(newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = list(reader.fieldnames)
        rows = list(reader)

    rows, fieldnames = update_runsheet(
        rows,
        fieldnames,
        mapping=mapping,
        glds=args.glds,
        entry_point=args.entry_point,
        mode=args.mode,
        assay_suffix=args.assay_suffix,
        strandedness=args.strandedness,
    )
    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
