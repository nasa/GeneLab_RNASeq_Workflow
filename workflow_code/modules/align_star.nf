process ALIGN_STAR {
  // Aligns reads against STAR index
  tag "Sample: ${ meta.id }"

  input:
    tuple val( meta ), path( reads )
    path(star_index_dir)

  output:
    path("${ meta.id }/${ meta.id }${ params.assay_suffix }*"), emit: publishables
    path("${ meta.id }/${ meta.id }${ params.assay_suffix }_Log.final.out"), emit: alignment_logs
    tuple val(meta), path("${ meta.id }/${ meta.id }${ params.assay_suffix }_Aligned.out.bam"), emit: bam
    tuple val(meta), path("${ meta.id }/${ meta.id }${ params.assay_suffix }_Aligned.toTranscriptome.out.bam"), emit: bam_to_transcriptome
    path("${ meta.id }/${ meta.id }${ params.assay_suffix }_ReadsPerGene.out.tab"), emit: reads_per_gene
    path("versions.yml"), emit: versions

  script:
    """
    STAR --twopassMode Basic \
    --outFilterType BySJout \
    --outSAMunmapped Within \
    --genomeDir ${ star_index_dir } \
    --outSAMattributes NH HI AS NM MD MC \
    --outFilterMismatchNoverReadLmax 0.04 \
    --outFilterMismatchNmax 999 \
    --outFilterMultimapNmax 20 \
    --alignIntronMin 20 \
    --alignSJoverhangMin 8 \
    ${ meta.paired_end ? "--alignMatesGapMax 1000000" : "" } \
    --alignIntronMax 1000000 \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 1 \
    --outSAMtype BAM Unsorted \
    --runThreadN ${ task.cpus } \
    --readFilesCommand zcat \
    --quantMode TranscriptomeSAM GeneCounts\
    --outFileNamePrefix '${ meta.id }/${ meta.id }${ params.assay_suffix }_' \
    --outReadsUnmapped Fastx \
    --genomeLoad NoSharedMemory \
    --readFilesIn ${ reads } 

    # Rename and gzip unmapped reads
    if [ ${meta.paired_end} == true ]; then
      mv "${ meta.id }/${ meta.id }${ params.assay_suffix }_Unmapped.out.mate1" "${ meta.id }/${ meta.id }${ params.assay_suffix }_R1_unmapped.fastq"
      mv "${ meta.id }/${ meta.id }${ params.assay_suffix }_Unmapped.out.mate2" "${ meta.id }/${ meta.id }${ params.assay_suffix }_R2_unmapped.fastq"
    else
      mv "${ meta.id }/${ meta.id }${ params.assay_suffix }_Unmapped.out.mate1" "${ meta.id }/${ meta.id }${ params.assay_suffix }_unmapped.fastq"
    fi
    gzip ${ meta.id }/${ meta.id }${ params.assay_suffix }*_unmapped.fastq

    echo '"${task.process}":' > versions.yml
    echo "    star: \$(STAR --version | sed 's/STAR_//')" >> versions.yml
    """
}
