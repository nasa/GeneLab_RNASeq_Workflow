/*
    Downloads BAM files for a specific sample ID from OSDR using file list TSV
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
    path(file_list)

    output:
    tuple val(meta), path("*.bam"), emit: bam_files

    script:
    def bam_suffix = params.mode == "microbes" ? 
        ".bam" :
        "_Aligned.toTranscriptome.out.bam"
    def output_name = params.mode == "microbes" ?
        "${meta.id}${params.assay_suffix}.bam" :
        "${meta.id}_Aligned.toTranscriptome.out.bam"
    
    """
    # Use provided file list TSV
    tsv_file="${file_list}"
    
    # Find BAM file for this sample - exact matching with assay suffix first
    bam_url=\$(grep "^${glds_accession}_rna_seq_${meta.id}${params.assay_suffix}${bam_suffix}" "\$tsv_file" | cut -f2 | head -1)
    
    # If not found, try without assay suffix (legacy naming)
    if [ -z "\$bam_url" ]; then
        bam_url=\$(grep "^${glds_accession}_rna_seq_${meta.id}${bam_suffix}" "\$tsv_file" | cut -f2 | head -1)
    fi
    
    if [ -z "\$bam_url" ]; then
        echo "ERROR: Could not find BAM file for ${meta.id} in file list"
        exit 1
    fi
    
    echo "Downloading BAM: \$bam_url"
    wget -q -O bam_temp.bam "\$bam_url" || exit 1
    
    mv bam_temp.bam "${output_name}"
    """
}
