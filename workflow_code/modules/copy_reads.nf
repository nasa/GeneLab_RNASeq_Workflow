process COPY_READS {
    tag "Sample: ${ meta.id }"

    publishDir "${ publishdir }/00-RawData/Fastq",
        pattern:  "*.gz" ,
        mode: params.publish_dir_mode
    // publishDir "${ publishdir }/00-RawData/Fastq",
    //     pattern:  "${meta.id}_R2_raw.fastq.gz" ,
    //     mode: params.publish_dir_mode
    // publishDir "${ publishdir }/00-RawData/Fastq",
    //     pattern:  "*.gz" ,
    //     mode: params.publish_dir_mode

    input:
        val(publishdir)
        tuple val(meta), path("?.gz")

    output:
        tuple val(meta), path("${meta.id}*.gz"), emit: raw_reads

    script:
        if ( meta.paired_end ) {
        """
        cp -P 1.gz ${meta.id}${params.assay_suffix}_R1_raw.fastq.gz
        cp -P 2.gz ${meta.id}${params.assay_suffix}_R2_raw.fastq.gz
        """
        } else {
        """
        cp -P 1.gz  ${meta.id}${params.assay_suffix}_raw.fastq.gz
        """
        }
}