/*
    Downloads read files for a specific sample ID from OSDR using the OSDR File Downloader
*/

process DOWNLOAD_OSDR_READS {
    tag "Sample: ${ meta.id }"
    
    errorStrategy 'ignore'  // Don't fail workflow if files not found
    
    publishDir { 
        read_type == "raw" ? 
            "${ publishdir }/00-RawData/Fastq" : 
            "${ publishdir }/01-TG_Preproc/Fastq" 
    },
        mode: params.publish_dir_mode

    input:
    val(publishdir)
    val(osd_accession)
    val(glds_accession)
    val(meta)
    val(read_type)  // "raw", "trimmed", etc.

    output:
    tuple val(meta), path("*.fastq.gz"), emit: trimmed_reads, optional: true
    path("*failed_download*.txt"), emit: failure_log, optional: true

    script:
    // Current naming with assay suffix
    def r1_current = "${glds_accession}_rna_seq_${meta.id}${params.assay_suffix}_R1_${read_type}.fastq.gz"
    def r2_current = "${glds_accession}_rna_seq_${meta.id}${params.assay_suffix}_R2_${read_type}.fastq.gz"
    def se_current = "${glds_accession}_rna_seq_${meta.id}${params.assay_suffix}_${read_type}.fastq.gz"
    // Legacy naming without assay suffix
    def r1_legacy = "${glds_accession}_rna_seq_${meta.id}_R1_${read_type}.fastq.gz"
    def r2_legacy = "${glds_accession}_rna_seq_${meta.id}_R2_${read_type}.fastq.gz"
    def se_legacy = "${glds_accession}_rna_seq_${meta.id}_${read_type}.fastq.gz"
    def se_legacy_r1 = "${glds_accession}_rna_seq_${meta.id}_R1_${read_type}.fastq.gz"
    
    """
    # Construct exact filenames based on GeneLab naming convention
    if [ "${meta.paired_end}" = "true" ]; then
        # Paired-end: try current naming first, then legacy
        echo "Trying current paired-end naming with assay suffix..."
        # Download R1 and R2 separately to avoid API issues with multiple search terms
        python3 ${projectDir}/bin/osdr_downloader.py \\
            --osd ${osd_accession} \\
            --measurement "transcription profiling" \\
            --tech "RNA-Seq" \\
            --search "${r1_current}" \\
            --out sample_downloads_current_r1
        
        python3 ${projectDir}/bin/osdr_downloader.py \\
            --osd ${osd_accession} \\
            --measurement "transcription profiling" \\
            --tech "RNA-Seq" \\
            --search "${r2_current}" \\
            --out sample_downloads_current_r2
        
        find sample_downloads_current_r1 -name "${r1_current}" -exec mv {} . \\; 2>/dev/null || true
        find sample_downloads_current_r2 -name "${r2_current}" -exec mv {} . \\; 2>/dev/null || true
        
        if [ ! -f "${r1_current}" ] || [ ! -f "${r2_current}" ]; then
            echo "Current PE files not found, trying legacy naming..."
            # Download R1 and R2 separately
            python3 ${projectDir}/bin/osdr_downloader.py \\
                --osd ${osd_accession} \\
                --measurement "transcription profiling" \\
                --tech "RNA-Seq" \\
                --search "${r1_legacy}" \\
                --out sample_downloads_legacy_r1
            
            python3 ${projectDir}/bin/osdr_downloader.py \\
                --osd ${osd_accession} \\
                --measurement "transcription profiling" \\
                --tech "RNA-Seq" \\
                --search "${r2_legacy}" \\
                --out sample_downloads_legacy_r2
            
            find sample_downloads_legacy_r1 -name "${r1_legacy}" -exec mv {} . \\; 2>/dev/null || true
            find sample_downloads_legacy_r2 -name "${r2_legacy}" -exec mv {} . \\; 2>/dev/null || true
            
            if [ ! -f "${r1_legacy}" ] || [ ! -f "${r2_legacy}" ]; then
                echo "WARNING: Could not find paired-end files for sample ${meta.id}"
                echo "Tried: ${r1_current}, ${r2_current}, ${r1_legacy}, ${r2_legacy}"
                echo "This dataset may not have trimmed reads available in OSDR."
                
                # Create failure log file
                cat > "${meta.id}${params.assay_suffix}_failed_download_${read_type}.txt" << EOF
Failed to download ${read_type} reads for sample: ${meta.id}
OSD: ${osd_accession}
GLDS: ${glds_accession}
Assay Suffix: ${params.assay_suffix}
Read Type: ${read_type}
Paired End: ${meta.paired_end}

Attempted file names:
- ${r1_current}
- ${r2_current}  
- ${r1_legacy}
- ${r2_legacy}

Reason: Files not found in OSDR repository
Date: \$(date)
EOF
                echo "Created failure log: ${meta.id}${params.assay_suffix}_failed_download_${read_type}.txt"
                exit 0  # Exit successfully to continue workflow
            fi
        fi
        
    else
        # Single-end: try current naming first, then legacy variants
        echo "Trying current single-end naming with assay suffix..."
        python3 ${projectDir}/bin/osdr_downloader.py \\
            --osd ${osd_accession} \\
            --measurement "transcription profiling" \\
            --tech "RNA-Seq" \\
            --search "${se_current}" \\
            --out sample_downloads_current
        
        find sample_downloads_current -name "${se_current}" -exec mv {} . \\; 2>/dev/null || true
        
        if [ ! -f "${se_current}" ]; then
            echo "Current SE file not found, trying legacy naming..."
            python3 ${projectDir}/bin/osdr_downloader.py \\
                --osd ${osd_accession} \\
                --measurement "transcription profiling" \\
                --tech "RNA-Seq" \\
                --search "${se_legacy}" \\
                --out sample_downloads_legacy
            
            find sample_downloads_legacy -name "${se_legacy}" -exec mv {} . \\; 2>/dev/null || true
            
            if [ ! -f "${se_legacy}" ]; then
                echo "Legacy SE file not found, trying legacy R1 naming..."
                python3 ${projectDir}/bin/osdr_downloader.py \\
                    --osd ${osd_accession} \\
                    --measurement "transcription profiling" \\
                    --tech "RNA-Seq" \\
                    --search "${se_legacy_r1}" \\
                    --out sample_downloads_legacy_r1
                
                find sample_downloads_legacy_r1 -name "${se_legacy_r1}" -exec mv {} . \\; 2>/dev/null || true
                
                if [ ! -f "${se_legacy_r1}" ]; then
                    echo "WARNING: Could not find single-end file for sample ${meta.id}"
                    echo "Tried: ${se_current}, ${se_legacy}, ${se_legacy_r1}"
                    echo "This dataset may not have trimmed reads available in OSDR."
                    
                    # Create failure log file
                    cat > "${meta.id}${params.assay_suffix}_failed_download_${read_type}.txt" << EOF
Failed to download ${read_type} reads for sample: ${meta.id}
OSD: ${osd_accession}
GLDS: ${glds_accession}
Assay Suffix: ${params.assay_suffix}
Read Type: ${read_type}
Paired End: ${meta.paired_end}

Attempted file names:
- ${se_current}
- ${se_legacy}
- ${se_legacy_r1}

Reason: Files not found in OSDR repository
Date: \$(date)
EOF
                    echo "Created failure log: ${meta.id}${params.assay_suffix}_failed_download_${read_type}.txt"
                    exit 0  # Exit successfully to continue workflow
                fi
            fi
        fi
    fi
    
    # Rename files to standardized names for downstream processing
    if [ "${meta.paired_end}" = "true" ]; then
        # Find the actual downloaded files and rename to standard format
        r1_file=\$(ls *R1*${read_type}.fastq.gz 2>/dev/null | head -1)
        r2_file=\$(ls *R2*${read_type}.fastq.gz 2>/dev/null | head -1)
        
        if [ -n "\$r1_file" ] && [ -n "\$r2_file" ]; then
            mv "\$r1_file" "${meta.id}_R1_${read_type}.fastq.gz"
            mv "\$r2_file" "${meta.id}_R2_${read_type}.fastq.gz"
            echo "Renamed paired-end files to standard format:"
            echo "  \$r1_file -> ${meta.id}_R1_${read_type}.fastq.gz"
            echo "  \$r2_file -> ${meta.id}_R2_${read_type}.fastq.gz"
        fi
    else
        # Single-end: rename to standard format
        se_file=\$(ls *${read_type}.fastq.gz 2>/dev/null | head -1)
        
        if [ -n "\$se_file" ]; then
            mv "\$se_file" "${meta.id}_${read_type}.fastq.gz"
            echo "Renamed single-end file to standard format:"
            echo "  \$se_file -> ${meta.id}_${read_type}.fastq.gz"
        fi
    fi
    
    echo "Final files for sample ${meta.id}:"
    ls -la *.fastq.gz
    """
}
