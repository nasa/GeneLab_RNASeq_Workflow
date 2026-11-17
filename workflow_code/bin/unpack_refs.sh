#!/usr/bin/env bash
set -euo pipefail
mkdir -p temp_fasta temp_gtf
cp -L "$1" temp_fasta/
cp -L "$2" temp_gtf/
shopt -s nullglob
if [[ -n $(echo temp_fasta/*.gz) ]]; then
    gunzip temp_fasta/*.gz
fi
if [[ -n $(echo temp_gtf/*.gz) ]]; then
    gunzip temp_gtf/*.gz
fi
for f in temp_fasta/*.fa temp_fasta/*.fna; do
    mv "$f" ./
done
for f in temp_gtf/*.gtf; do
    mv "$f" ./
done
rm -rf temp_fasta temp_gtf
