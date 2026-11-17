process VV_FASTQ {
  tag { fastq.name }
  label 'VV'

  input:
    path(fastq)

  output:
    path("*.vv_fastq.tsv"), emit: report

  when:
    !params.skip_vv

  script:
    """
    vv_check_fastq.sh "${fastq.toString()}" > "${fastq.name}.vv_fastq.tsv"
    """
}

process VV_RAW_READS {
  label 'VV'

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)
    path(raw_reads_multiqc_report_zip)
    path(fastq_reports, stageAs: 'fastq_reports/*')

  output:
    path("VV_log.csv"),   optional: params.skip_vv, emit: log
    path("versions.yml"), emit: versions

  script:
    def assay_suffix_arg = params.assay_suffix ? "--assay-suffix ${params.assay_suffix}" : ""
    def reports_arg = params.skip_vv ? "" : "--fastq-check-dir fastq_reports"
    """
    if ${ !params.skip_vv } ; then
      vv_raw_reads.py --runsheet ${runsheet} --outdir ${publishdir} ${assay_suffix_arg} ${reports_arg}
    fi

    echo '"${task.process}":' > versions.yml
    echo "    dp_tools: \$(pip show dp_tools | grep Version | sed 's/Version: //')" >> versions.yml
    """
}

process VV_TRIMMED_READS {
  label 'VV'

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)
    path(trimmed_reads_multiqc_report_zip)
    path(fastq_reports, stageAs: 'fastq_reports/*')

  output:
    path("VV_log.csv"), optional: params.skip_vv, emit: log
    path("versions.yml"), emit: versions

  script:
    def assay_suffix_arg = params.assay_suffix ? "--assay-suffix ${params.assay_suffix}" : ""
    def reports_arg = params.skip_vv ? "" : "--fastq-check-dir fastq_reports"
    """
    if ${ !params.skip_vv } ; then
      vv_trimmed_reads.py --runsheet ${runsheet} --outdir ${publishdir} ${assay_suffix_arg} ${reports_arg}
    fi

    echo '"${task.process}":' > versions.yml
    echo "    dp_tools: \$(pip show dp_tools | grep Version | sed 's/Version: //')" >> versions.yml
    """
}

process VV_BOWTIE2_ALIGNMENT {
  label 'VV'

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)  
    path(bowtie2_alignment_multiqc_report_zip)
    path(sorted_bams)
    
  output:
    path("VV_log.csv"), optional: params.skip_vv, emit: log
    path("versions.yml"), emit: versions

  script:
    def assay_suffix_arg = params.assay_suffix ? "--assay-suffix ${params.assay_suffix}" : ""
    """
    # Run V&V unless user requests to skip V&V
    if ${ !params.skip_vv } ; then
      vv_bowtie2_alignment.py --runsheet ${runsheet} --outdir ${publishdir} ${assay_suffix_arg}
    fi

    echo '"${task.process}":' > versions.yml
    echo "    dp_tools: \$(pip show dp_tools | grep Version | sed 's/Version: //')" >> versions.yml
    """
} 

process VV_RSEQC {
  label 'VV'

  input:
      path(dp_tools__NF_RCP)
      val(publishdir)
      val(meta)
      path(runsheet)
      path(rseqc_genebody_coverage_multiqc_report_zip)
      path(rseqc_infer_experiment_multiqc_report_zip)
      path(rseqc_inner_distance_multiqc_report_zip)
      path(rseqc_read_distribution_multiqc_report_zip)

  output:
      path("VV_log.csv"), optional: params.skip_vv, emit: log
      path("versions.yml"), emit: versions

  script:
    def assay_suffix_arg = params.assay_suffix ? "--assay_suffix ${params.assay_suffix}" : ""
    """
    # Run V&V unless user requests to skip V&V
    if ${ !params.skip_vv } ; then
      vv_rseqc.py --runsheet ${runsheet} --outdir ${publishdir} ${assay_suffix_arg}
    fi

    echo '"${task.process}":' > versions.yml
    echo "    dp_tools: \$(pip show dp_tools | grep Version | sed 's/Version: //')" >> versions.yml
    """
}

process VV_FEATURECOUNTS {
  label 'VV'

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)
    path(featurecounts_multiqc_report_zip)
    path(featurecounts_counts_rrnarm) // featurecounts counts (rRNArm)

  output:
    path("VV_log.csv"), optional: params.skip_vv, emit: log
    path("versions.yml"), emit: versions

  script:
  def assay_suffix_arg = params.assay_suffix ? "--assay-suffix ${params.assay_suffix}" : ""
  """
  # Run V&V unless user requests to skip V&V
  if ${ !params.skip_vv } ; then
    vv_featurecounts.py --runsheet ${runsheet} --outdir ${publishdir} ${assay_suffix_arg}
  fi

  echo '"${task.process}":' > versions.yml
  echo "    dp_tools: \$(pip show dp_tools | grep Version | sed 's/Version: //')" >> versions.yml
  """
}

process VV_DGE_DESEQ2 {
  label 'VV'

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)
    path(dge_table) // annotated dge table
    path(dge_table_rrnarm) // (rrna rm) annotated dge table
    val(stratify_by)

  output:
    path("VV_log.csv"), optional: params.skip_vv, emit: log
    path("log2fc_flag_characterization*.csv"), optional: true, emit: log2fc_flags
    path("versions.yml"), emit: versions
    
  script:
  def mode_params = params.mode == "microbes" ? "--mode microbes" : ""
  def assay_suffix_arg = params.assay_suffix ? "--assay_suffix ${params.assay_suffix}" : ""
  def stratify_arg = stratify_by ? "--stratify_by '${stratify_by}'" : ""
  """
  if ${ !params.skip_vv } ; then
    vv_dge_deseq2.py --runsheet ${runsheet} --outdir ${publishdir} ${assay_suffix_arg} ${mode_params} ${stratify_arg}
  fi

  echo '"${task.process}":' > versions.yml
  echo "    dp_tools: \$(pip show dp_tools | grep Version | sed 's/Version: //')" >> versions.yml
  """
}

process VV_STAR_ALIGNMENT {
  label 'VV'

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)
    path(star_alignment_multiqc_report_zip) 
    path(star_unnormalized_counts) 
    path(sorted_bams)
    
  output:
    path("VV_log.csv"), optional: params.skip_vv, emit: log

  script:
    def assay_suffix_arg = params.assay_suffix ? "--assay-suffix ${params.assay_suffix}" : ""
    """
    # Run V&V unless user requests to skip V&V
    if ${ !params.skip_vv } ; then
      vv_star_alignment.py --runsheet ${runsheet} --outdir ${publishdir} ${assay_suffix_arg}
    fi
    """
}

process VV_RSEM_COUNTS {
  label 'VV'

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)
    path(rsem_counts_multiqc_report_zip)
    path(rsem_genes_results_rrnarm) // RSEM sample.genes.results (rRNArm) 
    path(rsem_unnormalized_counts) // RSEM unnormalized counts
    
  output:
    path("VV_log.csv"), optional: params.skip_vv, emit: log
    path("versions.yml"), emit: versions
  
  script:
    def assay_suffix_arg = params.assay_suffix ? "--assay-suffix ${params.assay_suffix}" : ""
    """
    # Run V&V unless user requests to skip V&V
    if ${ !params.skip_vv } ; then
      vv_rsem_counts.py --runsheet ${runsheet} --outdir ${publishdir} ${assay_suffix_arg}
    fi

    echo '"${task.process}":' > versions.yml
    echo "    dp_tools: \$(pip show dp_tools | grep Version | sed 's/Version: //')" >> versions.yml
    """
}

process VV_CONCAT_FILTER {
  label 'VV'

  input:
    path("VV_in.csv")

  output:
    tuple path("VV_log_final${params.assay_suffix}.csv"), path("VV_log_final_only_issues${params.assay_suffix}.csv")

  script:
    def assay_suffix_arg = params.assay_suffix ? "--assay_suffix ${params.assay_suffix}" : ""
    """
    concat_logs.py ${assay_suffix_arg}
    filter_to_only_issues.py ${assay_suffix_arg}
    """
}