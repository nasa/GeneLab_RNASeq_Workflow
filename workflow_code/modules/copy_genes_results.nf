process COPY_GENES_RESULTS {
    tag "Sample: ${meta.id}"
    
    publishDir "${ publishdir }/03-RSEM_Counts/${meta.id}",
        mode: params.publish_dir_mode

    input:
        val(publishdir)
        tuple val(meta), path("input.genes.results")

    output:
        tuple val(meta), path("${meta.id}*.bam"), emit: genes_results

    script:
        def genes_results_name = "${meta.id}${params.assay_suffix}.genes.results"
        """
        # Copy and rename genes.results file to standard format
        cp "input.genes.results" "${genes_results_name}"
        echo "Copied genes.results: input.genes.results -> ${genes_results_name}"
        """
}

