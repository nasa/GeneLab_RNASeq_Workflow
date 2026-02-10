process COPY_GENES_RESULTS {
    tag "Sample: ${meta.id}"
    
    publishDir "${ publishdir }/03-RSEM_Counts/${meta.id}",
        mode: params.publish_dir_mode

    input:
        val(publishdir)
        tuple val(meta), path("?.genes.results")

    output:
        tuple val(meta), path("${meta.id}${params.assay_suffix}.genes.results"), emit: genes_results

    script:
        def genes_results_name = "${meta.id}${params.assay_suffix}.genes.results"
        """
        cp -P 1.genes.results ${genes_results_name}
        """
}

