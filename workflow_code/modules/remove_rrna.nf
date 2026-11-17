process REMOVE_RRNA {
    /**
     * Filters a RSEM `genes.results` file by removing rRNA gene IDs.
     */

    input:
        path rrna_ids_file
        tuple val(meta), path("${meta.id}${ params.assay_suffix }.genes.results")

    output:
        tuple val(meta), path("${meta.id}${params.assay_suffix}_rRNArm.genes.results"), emit: genes_results_rrnarm

    script:
        """
        sample_id="${meta.id}"
        filtered_file="\${sample_id}${params.assay_suffix}_rRNArm.genes.results"

        if [[ ! -s ${rrna_ids_file} ]]; then
            echo "WARNING: No rRNA gene IDs found for \${sample_id}; copying unfiltered genes.results to rRNArm output" >&2
            cp "${meta.id}${ params.assay_suffix }.genes.results" \${filtered_file}
        else
            awk 'NR==FNR {ids[\$1]=1; next} !(\$1 in ids)' ${rrna_ids_file} "${meta.id}${ params.assay_suffix }.genes.results" > \${filtered_file}
        fi
        """
}
