/*
    Downloads unnormalized counts table from OSDR using the OSDR File Downloader
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

    output:
    path("*Unnormalized_Counts${params.assay_suffix}.csv"), emit: counts_table, optional: true
    path("counts_table_failure${params.assay_suffix}.txt"), emit: failure_log, optional: true

    script:
    // Detect mode and set appropriate counts table filename
    def counts_current = params.mode == "microbes" ?
        "${glds_accession}_rna_seq_FeatureCounts_Unnormalized_Counts${params.assay_suffix}.csv" :
        "${glds_accession}_rna_seq_RSEM_Unnormalized_Counts${params.assay_suffix}.csv"
    def counts_legacy = params.mode == "microbes" ?
        "${glds_accession}_rna_seq_FeatureCounts_Unnormalized_Counts.csv" :
        "${glds_accession}_rna_seq_RSEM_Unnormalized_Counts.csv"
    
    """
    # Download unnormalized counts table
    echo "Trying current counts table naming with assay suffix..."
    python3 ${projectDir}/bin/osdr_downloader.py \\
        --osd ${osd_accession} \\
        --measurement "transcription profiling" \\
        --tech "RNA-Seq" \\
        --search "${counts_current}" \\
        --out sample_downloads_current
    
    find sample_downloads_current -name "${counts_current}" -exec mv {} . \\; 2>/dev/null || true
    
    if [ ! -f "${counts_current}" ]; then
        echo "Current counts table not found, trying without assay suffix..."
        python3 ${projectDir}/bin/osdr_downloader.py \\
            --osd ${osd_accession} \\
            --measurement "transcription profiling" \\
            --tech "RNA-Seq" \\
            --search "${counts_legacy}" \\
            --out sample_downloads_legacy
        
        find sample_downloads_legacy -name "${counts_legacy}" -exec mv {} . \\; 2>/dev/null || true
        
        if [ ! -f "${counts_legacy}" ]; then
            echo "WARNING: Could not find counts table"
            echo "Tried: ${counts_current}, ${counts_legacy}"
            echo "This dataset may not have counts table available in OSDR."
            
            # Create failure log file
            cat > "counts_table_failure${params.assay_suffix}.txt" << EOF
Failed to download counts table
OSD: ${osd_accession}
GLDS: ${glds_accession}
Assay Suffix: ${params.assay_suffix}
Mode: ${params.mode}

Attempted file names:
- ${counts_current}
- ${counts_legacy}

Reason: Files not found in OSDR
Date: \$(date)
EOF
            echo "Created failure log: counts_table_failure${params.assay_suffix}.txt"
        fi
    fi
    
    # Rename downloaded file to expected name
    if [ -f "${counts_current}" ]; then
        if [ "${params.mode}" == "microbes" ]; then
            mv "${counts_current}" "FeatureCounts_Unnormalized_Counts${params.assay_suffix}.csv"
        else
            mv "${counts_current}" "RSEM_Unnormalized_Counts${params.assay_suffix}.csv"
        fi
    elif [ -f "${counts_legacy}" ]; then
        if [ "${params.mode}" == "microbes" ]; then
            mv "${counts_legacy}" "FeatureCounts_Unnormalized_Counts${params.assay_suffix}.csv"
        else
            mv "${counts_legacy}" "RSEM_Unnormalized_Counts${params.assay_suffix}.csv"
        fi
    fi
    
    echo "Final files:"
    ls -la *.csv 2>/dev/null || echo "No counts table files found"
    """
}

