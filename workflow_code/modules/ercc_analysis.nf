process ERCC_ANALYSIS {
    input:
        val(has_ercc)
        val(accession)
        path(runsheet)
        val(assay_suffix)
        path(counts)
        val(organism)
        path("combined_ercc_analysis.ipynb")

    output:
        path("ercc_kernels.log")
        path("ercc_inputs.txt")
        path("combined_ercc_analysis${params.assay_suffix}.ipynb"), optional: true, emit: notebook
        path("ERCC_analysis${params.assay_suffix}.html"), optional: true, emit: html
        path("ERCC_analysis"), optional: true, emit: tables
        path("ERCC_analysis_error${params.assay_suffix}.txt"), optional: true, emit: error

    when:
        has_ercc

    script:
        def org = organism.toString().replaceAll(" ", "_").toLowerCase()
        def gene_id_prefix = 'EN'
        if ( org == 'arabidopsis_thaliana' ) {
            gene_id_prefix = 'AT'
        } else if ( org == 'rattus_norvegicus' ) {
            gene_id_prefix = 'ENSRN'
        } else if ( org == 'drosophila_melanogaster' ) {
            gene_id_prefix = 'FBgn'
        } else if ( org == 'homo_sapiens' ) {
            gene_id_prefix = 'ENSG'
        } else if ( org == 'brachypodium_distachyon' ) {
            gene_id_prefix = 'BRADI'
        }
        """
        export HOME="\$PWD"
        export JUPYTER_DATA_DIR="\$PWD/.jupyter"
        export JUPYTER_RUNTIME_DIR="\$PWD/.jupyter/runtime"
        mkdir -p "\$JUPYTER_RUNTIME_DIR"
        ercc_kernels.py \\
            --accession '${accession}' \\
            --runsheet '${runsheet}' \\
            --assay-suffix '${assay_suffix}' \\
            --counts '${counts}' \\
            --gene-id-prefix '${gene_id_prefix}' \\
            --notebook combined_ercc_analysis.ipynb \\
            | tee ercc_kernels.log
        """
}
