process UPDATE_OSDR_RUNSHEET {
    tag "${glds_accession}_${entry_point}"

    input:
        path(runsheet)
        path(file_list)
        val(glds_accession)
        val(entry_point)

    output:
        path("updated/*.csv"), emit: runsheet

    script:
        """
        mkdir -p updated
        update_osdr_runsheet.py \\
            --runsheet ${runsheet} \\
            --file-list ${file_list} \\
            --glds ${glds_accession} \\
            --entry-point ${entry_point} \\
            --strandedness ${params.strandedness} \\
            --mode ${params.mode} \\
            --assay-suffix '${params.assay_suffix}' \\
            --out 'updated/${runsheet.name}'
        """
}
