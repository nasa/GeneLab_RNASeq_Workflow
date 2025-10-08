/*
    Downloads dge table from OSDR using the OSDR File Downloader
*/

process DOWNLOAD_OSDR_DGE_TABLE {
    tag "Dataset: ${ glds_accession }"

    input:
    val(publishdir)
    val(osd_accession)
    val(glds_accession)

    output:
    path("*differential_expression*.csv"), emit: dge_table, optional: true
    path("dge_table_failure${params.assay_suffix}.txt"), emit: failure_log, optional: true

    script:
    // Set current (with suffix) and legacy (no suffix) DGE table filenames
    def dge_current = "${glds_accession}_rna_seq_differential_expression${params.assay_suffix}.csv"
    def dge_legacy  = "${glds_accession}_rna_seq_differential_expression.csv"
    
    """
    # Download DGE table
    echo "Trying current DGE table naming with assay suffix..."
    python3 ${projectDir}/bin/osdr_downloader.py \\
        --osd ${osd_accession} \\
        --measurement "transcription profiling" \\
        --tech "RNA-Seq" \\
        --search "${dge_current}" \\
        --out sample_downloads_current
    
    find sample_downloads_current -name "${dge_current}" -exec mv {} . \\; 2>/dev/null || true
    
    # Rename to expected name if found
    if [ -f "${dge_current}" ]; then
        mv "${dge_current}" "differential_expression${params.assay_suffix}.csv"
    fi
    
    if [ ! -f "differential_expression${params.assay_suffix}.csv" ]; then
        echo "Current DGE table not found, trying without assay suffix..."
        python3 ${projectDir}/bin/osdr_downloader.py \\
            --osd ${osd_accession} \\
            --measurement "transcription profiling" \\
            --tech "RNA-Seq" \\
            --search "${dge_legacy}" \\
            --out sample_downloads_legacy
        
        find sample_downloads_legacy -name "${dge_legacy}" -exec mv {} . \\; 2>/dev/null || true
        
        # Rename to expected name if found
        if [ -f "${dge_legacy}" ]; then
            mv "${dge_legacy}" "differential_expression${params.assay_suffix}.csv"
        fi
        
        if [ ! -f "differential_expression${params.assay_suffix}.csv" ]; then
            echo "WARNING: Could not find DGE table"
            echo "Tried: ${dge_current}, ${dge_legacy}"
            echo "This dataset may not have DGE table available in OSDR."
            
            # Create failure log file
            cat > "dge_table_failure${params.assay_suffix}.txt" << EOF
Failed to download DGE table
OSD: ${osd_accession}
GLDS: ${glds_accession}
Assay Suffix: ${params.assay_suffix}

Attempted file names:
- ${dge_current}
- ${dge_legacy}

Reason: Files not found in OSDR
Date: \$(date)
EOF
            echo "Created failure log: dge_table_failure${params.assay_suffix}.txt"
        fi
    fi
    
    echo "Final files:"
    ls -la *.csv 2>/dev/null || echo "No DGE table files found"
    """
}

