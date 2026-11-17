process TRIMGALORE {
  tag "Sample: ${ meta.id }"

  input:
    tuple val(meta), path(reads)

  output:
    tuple val(meta), path("${ meta.id }*trimmed.fastq.gz"), emit: reads
    path("${ meta.id }*.txt"), emit: reports
    path("versions.yml"), emit: versions

  script:
   // see https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md for --cores info
   def cores = (task.cpus && task.cpus <= 4) ? task.cpus : Math.min(task.cpus, 4)

    """
    trim_galore --gzip \
    --cores $cores \
    --phred33 \
    ${ meta.paired_end ? '--paired' : '' } \
    $reads \
    --output_dir .

    # rename with trimmed fastq files
    ${ meta.paired_end ? \
      "mv ${ meta.id }${ params.assay_suffix }_R1_raw_val_1.fq.gz ${ meta.id }${ params.assay_suffix }_R1_trimmed.fastq.gz; \
      mv ${ meta.id }${ params.assay_suffix }_R2_raw_val_2.fq.gz ${ meta.id }${ params.assay_suffix }_R2_trimmed.fastq.gz" : \
      "mv ${ meta.id }${ params.assay_suffix }_raw_trimmed.fq.gz ${ meta.id }${ params.assay_suffix }_trimmed.fastq.gz"}

    # rename trimming report files
    ${ meta.paired_end ? \
      "mv ${ meta.id }${ params.assay_suffix }_R1_raw.fastq.gz_trimming_report.txt ${ meta.id }_R1${ params.assay_suffix }_trimming_report.txt; \
      mv ${ meta.id }${ params.assay_suffix }_R2_raw.fastq.gz_trimming_report.txt ${ meta.id }_R2${ params.assay_suffix }_trimming_report.txt" : \
      "mv ${ meta.id }${ params.assay_suffix }_raw.fastq.gz_trimming_report.txt ${ meta.id }${ params.assay_suffix }_trimming_report.txt"}

    echo '"${task.process}":' > versions.yml
    echo "    trimgalore: \$(trim_galore -v | awk '/version/{print \$2}' | head -n1)" >> versions.yml
    echo "    cutadapt: \$(cutadapt --version)" >> versions.yml
    """
}
