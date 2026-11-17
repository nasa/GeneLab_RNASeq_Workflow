process FILTER_UNALIGNABLE {
    input:
        path(runsheet)
        path(mqc_data)

    output:
        path("dge_runsheet.csv"), emit: runsheet
        path("dropped_unalignable.tsv"), emit: dropped

    script:
        def mqc_arg = mqc_data.name != "NO_FILE" ? "--multiqc-data ${mqc_data}" : ""
        def assay_arg = params.assay_suffix ? "--assay-suffix ${params.assay_suffix}" : ""
        """
        filter_unalignable.py \\
          --runsheet ${runsheet} \\
          --threshold ${params.unalignable_threshold} \\
          ${assay_arg} \\
          ${mqc_arg}
        """
}
