/*
    Downloads DGE table from OSDR using file list TSV
*/

process DOWNLOAD_OSDR_DGE_TABLE {
    tag "Dataset: ${ glds_accession }"

    input:
    val(publishdir)
    val(osd_accession)
    val(glds_accession)
    path(file_list)

    output:
    path("*differential_expression*.csv"), emit: dge_table

    script:
    def output_name = "differential_expression${params.assay_suffix}.csv"
    
    """
    # Use provided file list TSV
    tsv_file="${file_list}"
    
    # Find DGE table - try with assay suffix first
    dge_url=\$(grep "^${glds_accession}_rna_seq_differential_expression${params.assay_suffix}.csv" "\$tsv_file" | cut -f2 | head -1)
    
    # If not found, try without assay suffix (legacy naming)
    if [ -z "\$dge_url" ]; then
        dge_url=\$(grep "^${glds_accession}_rna_seq_differential_expression.csv" "\$tsv_file" | cut -f2 | head -1)
    fi
    
    if [ -z "\$dge_url" ]; then
        echo "ERROR: Could not find differential expression table in file list"
        exit 1
    fi
    
    echo "Downloading DGE table: \$dge_url"
    wget -q -O dge_temp.csv "\$dge_url" || exit 1
    
    mv dge_temp.csv "${output_name}"
    """
}

