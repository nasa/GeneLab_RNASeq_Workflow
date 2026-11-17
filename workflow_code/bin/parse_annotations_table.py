#!/usr/bin/env python3
import pathlib
import shlex
import sys

csv_path = pathlib.Path(sys.argv[1])
organism_sci = sys.argv[2]
entry_point = sys.argv[3]
have_refs = sys.argv[4] == "1"

organisms = {}
for line in csv_path.read_text().splitlines():
    fields = line.split(",")
    if len(fields) > 1:
        organisms[fields[1]] = fields

key = (organism_sci[:1].upper() + organism_sci[1:]).replace("_", " ")


def write_env(fa, gtf, ga, src, ver):
    pathlib.Path("parsed.env").write_text(
        f"fasta_url={shlex.quote(fa)}\n"
        f"gtf_url={shlex.quote(gtf)}\n"
        f"gene_annotations_url={shlex.quote(ga)}\n"
        f"reference_source={shlex.quote(src)}\n"
        f"reference_version={shlex.quote(ver)}\n"
    )


def col(row, i):
    return row[i] if len(row) > i else ""


if key in organisms:
    row = organisms[key]
    write_env(col(row, 5), col(row, 6), col(row, 10), col(row, 4), col(row, 3))
    print(f"Reading annotations table from {csv_path.resolve()}")
    print(f"Annotation table values parsed for '{key}':")
    print(f"            Reference Fasta URL: {col(row, 5)}")
    print(f"            Reference GTF URL: {col(row, 6)}")
    print(f"            Gene Annotations URL: {col(row, 10)}")
    print(f"            Reference Source: {col(row, 4)}")
    if "ensembl" in col(row, 4).lower():
        print(f"            Reference Version: {col(row, 3)}")
    sys.exit(0)

if entry_point != "dge_table" and not have_refs:
    sys.exit(
        f"Organism '{key}' is not in the annotations table. "
        "Pass --reference_fasta and --reference_gtf."
    )

print(f"WARNING: Organism '{key}' not found in annotations table.")
print("Returning null values for all outputs.")
write_env("", "", "", "", "")
