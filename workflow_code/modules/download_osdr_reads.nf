/*
    Downloads read files for a specific sample ID from OSDR using file list TSV
*/

process DOWNLOAD_OSDR_READS {
    tag "Sample: ${ meta.id }"
    
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
    path(file_list)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: trimmed_reads

    script:
    """
    tsv_file="${file_list}"
    
    # Parse file list to find exact files for this sample and read type
    if [ "${meta.paired_end}" = "true" ]; then
        # Find R1 and R2 files - exact matching with assay suffix first
        r1_url=\$(grep "^${glds_accession}_rna.seq_${meta.id}${params.assay_suffix}_R1_${read_type}.fastq.gz" "\$tsv_file" | cut -f2 | head -1)
        r2_url=\$(grep "^${glds_accession}_rna.seq_${meta.id}${params.assay_suffix}_R2_${read_type}.fastq.gz" "\$tsv_file" | cut -f2 | head -1)
        
        # If not found, try without assay suffix (legacy naming)
        if [ -z "\$r1_url" ] || [ -z "\$r2_url" ]; then
            r1_url=\$(grep "^${glds_accession}_rna.seq_${meta.id}_R1_${read_type}.fastq.gz" "\$tsv_file" | cut -f2 | head -1)
            r2_url=\$(grep "^${glds_accession}_rna.seq_${meta.id}_R2_${read_type}.fastq.gz" "\$tsv_file" | cut -f2 | head -1)
        fi
        
        if [ -z "\$r1_url" ] || [ -z "\$r2_url" ]; then
            echo "ERROR: Could not find paired-end ${read_type} reads for ${meta.id} in file list"
            exit 1
        fi
        
        echo "Downloading R1: \$r1_url"
        wget -q -O r1.fastq.gz "\$r1_url" || exit 1
        
        echo "Downloading R2: \$r2_url"
        wget -q -O r2.fastq.gz "\$r2_url" || exit 1
        
        mv r1.fastq.gz "${meta.id}${params.assay_suffix}_R1_${read_type}.fastq.gz"
        mv r2.fastq.gz "${meta.id}${params.assay_suffix}_R2_${read_type}.fastq.gz"
    else
        # Find SE file - exact matching with assay suffix first
        se_url=\$(grep "^${glds_accession}_rna.seq_${meta.id}${params.assay_suffix}_${read_type}.fastq.gz" "\$tsv_file" | cut -f2 | head -1)
        
        # If not found, try with R1 suffix (some SE data has R1 in the name)
        if [ -z "\$se_url" ]; then
            se_url=\$(grep "^${glds_accession}_rna.seq_${meta.id}${params.assay_suffix}_R1_${read_type}.fastq.gz" "\$tsv_file" | cut -f2 | head -1)
        fi
        
        # If still not found, try without assay suffix (legacy naming)
        if [ -z "\$se_url" ]; then
            se_url=\$(grep "^${glds_accession}_rna.seq_${meta.id}_${read_type}.fastq.gz" "\$tsv_file" | cut -f2 | head -1)
        fi
        
        # If still not found, try legacy with R1
        if [ -z "\$se_url" ]; then
            se_url=\$(grep "^${glds_accession}_rna.seq_${meta.id}_R1_${read_type}.fastq.gz" "\$tsv_file" | cut -f2 | head -1)
        fi
        
        if [ -z "\$se_url" ]; then
            echo "ERROR: Could not find single-end ${read_type} reads for ${meta.id} in file list"
            exit 1
        fi
        
        echo "Downloading SE: \$se_url"
        wget -q -O se.fastq.gz "\$se_url" || exit 1
        
        mv se.fastq.gz "${meta.id}${params.assay_suffix}_${read_type}.fastq.gz"
    fi
    """
}
