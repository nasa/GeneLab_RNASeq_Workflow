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
    path("*Unnormalized_Counts*.csv"), emit: counts_table, optional: true
    path("*failed_download*.txt"), emit: failure_log, optional: true

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
        echo "Current counts table not found, trying legacy naming..."
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
            cat > "failed_download_counts_table.txt" << EOF
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
            echo "Created failure log: failed_download_counts_table.txt"
        fi
    fi
    
    # Rename to standard format (remove GLDS prefix)
    counts_file=\$(ls *Unnormalized_Counts*.csv 2>/dev/null | head -1)
    
    if [ -n "\$counts_file" ]; then
        standard_name=\$(echo "\$counts_file" | sed "s/${glds_accession}_rna_seq_//")
        mv "\$counts_file" "\$standard_name"
        echo "Renamed counts table to standard format:"
        echo "  \$counts_file -> \$standard_name"
    fi
    
    echo "Final files:"
    ls -la *.csv
    """
}

