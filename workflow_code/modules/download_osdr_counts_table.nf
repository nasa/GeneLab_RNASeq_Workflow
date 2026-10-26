/*
    Downloads unnormalized counts table from OSDR using file list TSV
*/

process DOWNLOAD_OSDR_COUNTS_TABLE {
    tag "Dataset: ${ glds_accession }"
    
    publishDir { 
        params.mode == "microbes" ? 
            "${ publishdir }/03-FeatureCounts" : 
            "${ publishdir }/03-RSEM_Counts" 
    },
        mode: params.publish_dir_mode

    input:
    val(publishdir)
    val(osd_accession)
    val(glds_accession)
    path(file_list)

    output:
    path("*Unnormalized_Counts${params.assay_suffix}.csv"), emit: counts_table

    script:
    def counts_type = params.mode == "microbes" ? "FeatureCounts" : "RSEM"
    def output_name = "${counts_type}_Unnormalized_Counts${params.assay_suffix}.csv"
    
    """
    # Use provided file list TSV
    tsv_file="${file_list}"
    
    # Find counts table - try with assay suffix first
    counts_url=\$(grep "^${glds_accession}_rna_seq_${counts_type}_Unnormalized_Counts${params.assay_suffix}.csv" "\$tsv_file" | cut -f2 | head -1)
    
    # If not found, try without assay suffix (legacy naming)
    if [ -z "\$counts_url" ]; then
        counts_url=\$(grep "^${glds_accession}_rna_seq_${counts_type}_Unnormalized_Counts.csv" "\$tsv_file" | cut -f2 | head -1)
    fi
    
    if [ -z "\$counts_url" ]; then
        echo "ERROR: Could not find ${counts_type} raw counts table in file list"
        exit 1
    fi
    
    echo "Downloading counts table: \$counts_url"
    wget -q -O counts_temp.csv "\$counts_url" || exit 1
    
    mv counts_temp.csv "${output_name}"
    """
}

