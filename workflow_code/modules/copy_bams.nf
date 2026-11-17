process COPY_BAMS {
    tag "Sample: ${meta.id}"

    input:
        tuple val(meta), path("?.bam")

    output:
        tuple val(meta), path("${meta.id}*.bam"), emit: bam_files

    script:
        def bam_name = meta.genome_bam ? "${meta.id}${params.assay_suffix}_Aligned.out.bam" : params.mode == "microbes" ?
            "${meta.id}${params.assay_suffix}.bam" :
            "${meta.id}${params.assay_suffix}_Aligned.toTranscriptome.out.bam"
        """
        cp -P 1.bam ${bam_name}
        """
}

