process COPY_READS {
    tag "Sample: ${ meta.id }"

    input:
        tuple val(meta), path("?.gz")
        val(type) // "raw", "trimmed"

    output:
        tuple val(meta), path("${meta.id}*.gz"), emit: reads

    script:
        def suffix = type == "raw" ? "raw" : "trimmed"
        def nrec = params.truncate_to ? params.truncate_to.toInteger() : 0
        def out1 = "${meta.id}${params.assay_suffix}_R1_${suffix}.fastq.gz"
        def out2 = "${meta.id}${params.assay_suffix}_R2_${suffix}.fastq.gz"
        def out_se = "${meta.id}${params.assay_suffix}_${suffix}.fastq.gz"
        if ( meta.paired_end && nrec ) {
        """
        fastq_first_n.sh ${nrec} ${out1} 1.gz
        fastq_first_n.sh ${nrec} ${out2} 2.gz
        """
        } else if ( meta.paired_end ) {
        """
        cp -P 1.gz ${out1}
        cp -P 2.gz ${out2}
        """
        } else if ( nrec ) {
        """
        fastq_first_n.sh ${nrec} ${out_se} 1.gz
        """
        } else {
        """
        cp -P 1.gz ${out_se}
        """
        }
}
