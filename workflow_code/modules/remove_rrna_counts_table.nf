process REMOVE_RRNA_COUNTS_TABLE {

  input:
    path(counts_table)
    path(rrna_ids)

  output:
    path("*Unnormalized_Counts_rRNArm*.csv"), emit: counts_rrnarm

  script:
    def output_name = params.mode == "microbes" ?
        "FeatureCounts_Unnormalized_Counts_rRNArm${params.assay_suffix}.csv" :
        "RSEM_Unnormalized_Counts_rRNArm${params.assay_suffix}.csv"
    """
    # Remove rRNA genes from counts table
    remove_rrna_counts_table.py \\
        --counts_table ${counts_table} \\
        --rrna_ids ${rrna_ids} \\
        --output ${output_name}
    """
}
