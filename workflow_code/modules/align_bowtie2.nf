/*
    Aligns reads against a Bowtie2 index and writes the SAM.
    SAM_TO_BAM converts it and drops unmapped reads.
*/

process ALIGN_BOWTIE2 {
  tag "Sample: ${ meta.id }"

  input:
    tuple val(meta), path(reads)
    path(bowtie2_index_dir)

  output:
    tuple val(meta), path("${ meta.id }${ params.assay_suffix }.sam"), emit: sam
    tuple val(meta), path("${ meta.id }${ params.assay_suffix }_*unmapped.fastq.gz"), emit: unmapped_reads
    tuple val(meta), path("${ meta.id }${ params.assay_suffix }.bowtie2.log"), emit: alignment_logs
    path("versions.yml"), emit: versions

  script:
    def unaligned = meta.paired_end ? "--un-conc-gz ${ meta.id }${ params.assay_suffix }.unmapped.fastq.gz" : "--un-gz ${ meta.id }${ params.assay_suffix }.unmapped.fastq.gz"
    def readArgs = meta.paired_end ? "-1 ${reads[0]} -2 ${reads[1]}" : "-U ${reads}"

    """
    INDEX=\$(find -L ${ bowtie2_index_dir } -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//")

    bowtie2 -x "\$INDEX" \\
      -a \\
      ${readArgs} \\
      --threads ${task.cpus} \\
      --minins 10 \\
      --maxins 1000 \\
      ${ unaligned } \\
      2> ${ meta.id }${ params.assay_suffix }.bowtie2.log \\
      > ${ meta.id }${ params.assay_suffix }.sam

    # Rename unmapped reads
    if [ ${ meta.paired_end } == true ]; then
      # For paired-end data
      mv "${ meta.id }${ params.assay_suffix }.unmapped.fastq.1.gz" "${ meta.id }${ params.assay_suffix }_R1_unmapped.fastq.gz"
      mv "${ meta.id }${ params.assay_suffix }.unmapped.fastq.2.gz" "${ meta.id }${ params.assay_suffix }_R2_unmapped.fastq.gz"
    else
      # For single-end data
      mv "${ meta.id }${ params.assay_suffix }.unmapped.fastq.gz" "${ meta.id }${ params.assay_suffix }_unmapped.fastq.gz"
    fi

    echo '"${ task.process }":' > versions.yml
    echo "    bowtie2: \$(bowtie2 --version | head -n1 | awk '{print \$3}')" >> versions.yml
    """
}
