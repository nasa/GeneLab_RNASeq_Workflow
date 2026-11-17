#!/usr/bin/env python
import sys
import pandas as pd

def main(counts_file, rrna_ids_file, filtered_output):
    with open(rrna_ids_file) as f:
        rrna_ids = set(line.strip() for line in f if line.strip())

    header_lines = []
    with open(counts_file) as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                break

    df = pd.read_csv(counts_file, sep='\t', comment='#')
    gene_id_col = df.columns[0]
    rna_mask = df[gene_id_col].isin(rrna_ids)

    with open(filtered_output, 'w') as f:
        for line in header_lines:
            f.write(line)
        df[~rna_mask].to_csv(f, sep='\t', index=False)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit("Usage: python remove_rrna_featurecounts.py <counts_file> <rrna_ids_file> <filtered_output>")
    main(sys.argv[1], sys.argv[2], sys.argv[3])
