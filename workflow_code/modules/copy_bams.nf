process COPY_BAMS {
    tag "Sample: ${meta.id}"
    
    publishDir { 
        params.mode == "microbes" ? 
            "${ publishdir }/02-Bowtie2_Alignment/${meta.id}" : 
            "${ publishdir }/02-STAR_Alignment/${meta.id}" 
    },
        mode: params.publish_dir_mode

    input:
        val(publishdir)
        tuple val(meta), path("?.bam")

    output:
        tuple val(meta), path("${meta.id}*.bam"), emit: bam_files

    script:
        def bam_name = params.mode == "microbes" ? 
            "${meta.id}${params.assay_suffix}.bam" :
            "${meta.id}${params.assay_suffix}_Aligned.toTranscriptome.out.bam"
        """
        # Copy and rename BAM file to standard format
        input_bam=\$(ls *.bam | head -1)
        cp "\$input_bam" "${bam_name}"
        echo "Copied BAM: \$input_bam -> ${bam_name}"
        """
}

