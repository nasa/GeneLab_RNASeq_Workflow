process COPY_COUNTS_TABLE {

    input:
        path("?.csv")

    output:
        path("*Unnormalized_Counts*.csv"), emit: counts_table

    script:
        def counts_name = params.mode == "microbes" ? 
            "FeatureCounts_Unnormalized_Counts${params.assay_suffix}.csv" :
            "RSEM_Unnormalized_Counts${params.assay_suffix}.csv"
        """
        cp -P 1.csv ${counts_name}
        """
}

