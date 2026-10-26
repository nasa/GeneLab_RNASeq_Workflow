process COPY_READS {
    tag "Sample: ${ meta.id }"

    publishDir "${ publishdir }/${ type == "raw" ? "00-RawData/Fastq" : "01-TG_Preproc/Fastq" }",
        pattern:  "*.gz" ,
        mode: params.publish_dir_mode

    input:
        val(publishdir)
        tuple val(meta), path("?.gz")
        val(type) // "raw", "trimmed"

    output:
        tuple val(meta), path("${meta.id}*.gz"), emit: reads

    script:
        def suffix = type == "raw" ? "raw" : "trimmed"
        if ( meta.paired_end ) {
        """
        cp -P 1.gz ${meta.id}${params.assay_suffix}_R1_${suffix}.fastq.gz
        cp -P 2.gz ${meta.id}${params.assay_suffix}_R2_${suffix}.fastq.gz
        """
        } else {
        """
        cp -P 1.gz  ${meta.id}${params.assay_suffix}_${suffix}.fastq.gz
        """
        }
}
