/*
    Downloads RSEM genes.results files for a specific sample ID from OSDR using the OSDR File Downloader
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

    output:
    tuple val(meta), path("*.genes.results"), emit: genes_results, optional: true
    path("*failed_download*.txt"), emit: failure_log, optional: true

    script:
    // Current naming with assay suffix
    def genes_results_current = "${glds_accession}_rna_seq_${meta.id}${params.assay_suffix}.genes.results"
    def genes_results_legacy = "${glds_accession}_rna_seq_${meta.id}.genes.results"
    
    """
    # Download RSEM genes.results file
    echo "Trying current genes.results naming with assay suffix..."
    python3 ${projectDir}/bin/osdr_downloader.py \\
        --osd ${osd_accession} \\
        --measurement "transcription profiling" \\
        --tech "RNA-Seq" \\
        --search "${genes_results_current}" \\
        --out sample_downloads_current
    
    find sample_downloads_current -name "${genes_results_current}" -exec mv {} . \\; 2>/dev/null || true
    
    if [ ! -f "${genes_results_current}" ]; then
        echo "Current genes.results file not found, trying legacy naming..."
        python3 ${projectDir}/bin/osdr_downloader.py \\
            --osd ${osd_accession} \\
            --measurement "transcription profiling" \\
            --tech "RNA-Seq" \\
            --search "${genes_results_legacy}" \\
            --out sample_downloads_legacy
        
        find sample_downloads_legacy -name "${genes_results_legacy}" -exec mv {} . \\; 2>/dev/null || true
        
        if [ ! -f "${genes_results_legacy}" ]; then
            echo "WARNING: Could not find genes.results file for sample ${meta.id}"
            echo "Tried: ${genes_results_current}, ${genes_results_legacy}"
            echo "This dataset may not have genes.results files available in OSDR."
            
            # Create failure log file
            cat > "${meta.id}${params.assay_suffix}_failed_download_genes_results.txt" << EOF
Failed to download genes.results file for sample: ${meta.id}
OSD: ${osd_accession}
GLDS: ${glds_accession}
Assay Suffix: ${params.assay_suffix}

Attempted file names:
- ${genes_results_current}
- ${genes_results_legacy}

Reason: Files not found in OSDR
Date: \$(date)
EOF
            echo "Created failure log: ${meta.id}${params.assay_suffix}_failed_download_genes_results.txt"
        fi
    fi
    
    # Rename to standardized name for downstream processing
    genes_file=\$(ls *.genes.results 2>/dev/null | head -1)
    
    if [ -n "\$genes_file" ]; then
        mv "\$genes_file" "${meta.id}${params.assay_suffix}.genes.results"
        echo "Renamed genes.results file to standard format:"
        echo "  \$genes_file -> ${meta.id}${params.assay_suffix}.genes.results"
    fi
    
    echo "Final files for sample ${meta.id}:"
    ls -la *.genes.results
    """
}