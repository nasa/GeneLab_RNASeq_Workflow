/*
    Downloads BAM files for a specific sample ID from OSDR using the OSDR File Downloader
*/

process DOWNLOAD_OSDR_BAM {
    tag "Sample: ${ meta.id }"
    
    publishDir "${ publishdir }/02-STAR_Alignment/${meta.id}",
        mode: params.publish_dir_mode

    input:
    val(publishdir)
    val(osd_accession)
    val(glds_accession)
    val(meta)

    output:
    tuple val(meta), path("*.bam"), emit: bam_files, optional: true
    path("*failed_download*.txt"), emit: failure_log, optional: true

    script:
    // Detect mode and set appropriate BAM filename pattern
    def bam_current = params.mode == "microbes" ? 
        "${glds_accession}_rna_seq_${meta.id}${params.assay_suffix}.bam" :
        "${glds_accession}_rna_seq_${meta.id}${params.assay_suffix}_Aligned.toTranscriptome.out.bam"
    def bam_legacy = params.mode == "microbes" ?
        "${glds_accession}_rna_seq_${meta.id}.bam" :
        "${glds_accession}_rna_seq_${meta.id}_Aligned.toTranscriptome.out.bam"
    
    """
    # Download BAM file for quantification
    echo "Trying current BAM naming with assay suffix..."
    python3 ${projectDir}/bin/osdr_downloader.py \\
        --osd ${osd_accession} \\
        --measurement "transcription profiling" \\
        --tech "RNA-Seq" \\
        --search "${bam_current}" \\
        --out sample_downloads_current
    
    find sample_downloads_current -name "${bam_current}" -exec mv {} . \\; 2>/dev/null || true
    
    if [ ! -f "${bam_current}" ]; then
        echo "Current BAM file not found, trying legacy naming..."
        python3 ${projectDir}/bin/osdr_downloader.py \\
            --osd ${osd_accession} \\
            --measurement "transcription profiling" \\
            --tech "RNA-Seq" \\
            --search "${bam_legacy}" \\
            --out sample_downloads_legacy
        
        find sample_downloads_legacy -name "${bam_legacy}" -exec mv {} . \\; 2>/dev/null || true
        
        if [ ! -f "${bam_legacy}" ]; then
            echo "WARNING: Could not find BAM file for sample ${meta.id}"
            echo "Tried: ${bam_current}, ${bam_legacy}"
            echo "This dataset may not have BAM files available in OSDR."
            
            # Create failure log file
            cat > "${meta.id}${params.assay_suffix}_failed_download_bam.txt" << EOF
Failed to download BAM file for sample: ${meta.id}
OSD: ${osd_accession}
GLDS: ${glds_accession}
Assay Suffix: ${params.assay_suffix}

Attempted file names:
- ${bam_current}
- ${bam_legacy}

Reason: Files not found in OSDR
Date: \$(date)
EOF
            echo "Created failure log: ${meta.id}${params.assay_suffix}_failed_download_bam.txt"
        fi
    fi
    
    # Rename to standardized name for downstream processing
    bam_file=\$(ls *.bam 2>/dev/null | head -1)
    
    if [ -n "\$bam_file" ]; then
        if [[ "${params.mode}" == "microbes" ]]; then
            mv "\$bam_file" "${meta.id}${params.assay_suffix}.bam"
            echo "Renamed BAM file to standard format:"
            echo "  \$bam_file -> ${meta.id}${params.assay_suffix}.bam"
        else
            mv "\$bam_file" "${meta.id}_Aligned.toTranscriptome.out.bam"
            echo "Renamed BAM file to standard format:"
            echo "  \$bam_file -> ${meta.id}_Aligned.toTranscriptome.out.bam"
        fi
    fi
    
    echo "Final files for sample ${meta.id}:"
    ls -la *.bam
    """
}
