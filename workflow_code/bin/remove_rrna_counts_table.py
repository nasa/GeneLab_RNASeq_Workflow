#!/usr/bin/env python 
"""
Remove a column list of genes (rRNA genes) from a counts table 
"""

import argparse
import pandas as pd
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description='Remove rRNA genes from counts table')
    parser.add_argument('--counts_table', required=True, help='Path to counts table CSV')
    parser.add_argument('--rrna_ids', required=True, help='Path to rRNA gene IDs file')
    parser.add_argument('--output', required=True, help='Output path for rRNA-removed counts table')
    
    args = parser.parse_args()
    
    # Read rRNA gene IDs
    rrna_genes = set()
    with open(args.rrna_ids, 'r') as f:
        for line in f:
            gene_id = line.strip()
            if gene_id:
                rrna_genes.add(gene_id)
    
    print(f"Loaded {len(rrna_genes)} rRNA gene IDs")
    if rrna_genes:
        print(f"Example rRNA IDs: {list(rrna_genes)[:5]}")
    
    # Read counts table
    counts_df = pd.read_csv(args.counts_table, index_col=0)
    original_genes = len(counts_df)
    print(f"Original counts table: {original_genes} genes")
    
    # Remove rRNA genes
    genes_to_remove = [gene for gene in counts_df.index if gene in rrna_genes]
    print(f"Removing {len(genes_to_remove)} rRNA genes")
    
    counts_df_filtered = counts_df[~counts_df.index.isin(rrna_genes)]
    remaining_genes = len(counts_df_filtered)
    print(f"Filtered counts table: {remaining_genes} genes")
    
    # Save filtered table
    counts_df_filtered.to_csv(args.output)
    print(f"Saved rRNA-removed counts table to: {args.output}")


if __name__ == "__main__":
    main()

