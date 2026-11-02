process COPY_COUNTS_TABLE {

    publishDir { 
        params.mode == "microbes" ? 
            "${ publishdir }/03-FeatureCounts" : 
            "${ publishdir }/03-RSEM_Counts" 
    },
        mode: params.publish_dir_mode

    input:
        val(publishdir)
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

