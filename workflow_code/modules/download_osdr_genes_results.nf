/*
    Downloads RSEM genes.results files for a specific sample ID from OSDR using file list TSV
*/

process DOWNLOAD_OSDR_GENES_RESULTS {
    tag "Sample: ${ meta.id }"
    
    publishDir "${ publishdir }/03-RSEM_Counts/${meta.id}",
        mode: params.publish_dir_mode

    input:
    val(publishdir)
    val(osd_accession)
    val(glds_accession)
    val(meta)
    path(file_list)

    output:
    tuple val(meta), path("*.genes.results"), emit: genes_results

    script:
    """
    # Use provided file list TSV
    tsv_file="${file_list}"
    
    # Find genes.results file for this sample - match exact filename
    genes_url=\$(grep "^${glds_accession}_rna_seq_${meta.id}${params.assay_suffix}.genes.results" "\$tsv_file" | cut -f2 | head -1)
    
    # If not found, try without assay suffix (legacy naming)
    if [ -z "\$genes_url" ]; then
        genes_url=\$(grep "^${glds_accession}_rna_seq_${meta.id}.genes.results" "\$tsv_file" | cut -f2 | head -1)
    fi
    
    if [ -z "\$genes_url" ]; then
        echo "ERROR: Could not find genes.results file for ${meta.id} in file list"
        exit 1
    fi
    
    echo "Downloading genes.results: \$genes_url"
    wget -q -O genes_temp.genes.results "\$genes_url" || exit 1
    
    mv genes_temp.genes.results "${meta.id}${params.assay_suffix}.genes.results"
    """
}