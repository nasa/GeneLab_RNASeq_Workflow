process VV_RAW_READS {
  // Run vv_raw_reads.py to generate VV_log.csv pointing to already published data
  label 'VV'

  // Publish VV log
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.tokenize(':').last() }${ params.assay_suffix }.csv" }

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)                
    path(raw_reads_multiqc_report_zip)

  output:
    path("VV_log.csv"),   optional: params.skip_vv, emit: log
    path("versions.yml"), emit: versions

  script:
    """
    # Run V&V unless user requests to skip V&V
    if ${ !params.skip_vv } ; then
      vv_raw_reads.py --runsheet ${runsheet} --outdir ${publishdir} --assay-suffix ${params.assay_suffix}
    fi

    echo '"${task.process}":' > versions.yml
    echo "    dp_tools: \$(pip show dp_tools | grep Version | sed 's/Version: //')" >> versions.yml
    """
}

process VV_TRIMMED_READS {
  label 'VV'

  // Log publishing
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.tokenize(':').last() }${ params.assay_suffix }.csv" }

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)                
    path(trimmed_reads_multiqc_report_zip)

  output:
    path("VV_log.csv"), optional: params.skip_vv, emit: log
    path("versions.yml"), emit: versions

  script:
    """
    # Run V&V unless user requests to skip V&V
    if ${ !params.skip_vv } ; then
      vv_trimmed_reads.py --runsheet ${runsheet} --outdir ${publishdir} --assay-suffix ${params.assay_suffix}
    fi

    echo '"${task.process}":' > versions.yml
    echo "    dp_tools: \$(pip show dp_tools | grep Version | sed 's/Version: //')" >> versions.yml
    """
}

process VV_BOWTIE2_ALIGNMENT {
  // Log publishing
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_VV_ALIGNMENT${ params.assay_suffix }.csv" }

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
    """
    # Run V&V unless user requests to skip V&V
    if ${ !params.skip_vv } ; then
      vv_bowtie2_alignment.py --runsheet ${runsheet} --outdir ${publishdir} --assay-suffix ${params.assay_suffix}
    fi

    echo '"${task.process}":' > versions.yml
    echo "    dp_tools: \$(pip show dp_tools | grep Version | sed 's/Version: //')" >> versions.yml
    """
} 

process VV_RSEQC {
  // Log publishing
  publishDir "${ publishdir }",
    pattern: "VV_log.csv",
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.tokenize(':').last() }${ params.assay_suffix }.csv" }

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
    """
    # Run V&V unless user requests to skip V&V
    if ${ !params.skip_vv } ; then
      vv_rseqc.py --runsheet ${runsheet} --outdir ${publishdir} --assay_suffix ${params.assay_suffix}
    fi

    echo '"${task.process}":' > versions.yml
    echo "    dp_tools: \$(pip show dp_tools | grep Version | sed 's/Version: //')" >> versions.yml
    """
}

process VV_FEATURECOUNTS {
  // Log publishing
  publishDir "${ publishdir }",
    pattern: "VV_log.csv",
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_VV_COUNTS${ params.assay_suffix }.csv" }

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

  script:
  """
  # Run V&V unless user requests to skip V&V
  if ${ !params.skip_vv } ; then
    vv_featurecounts.py --runsheet ${runsheet} --outdir ${publishdir} --assay-suffix ${params.assay_suffix}
  fi
  """
}

process VV_DGE_DESEQ2 {
  // Log publishing
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_VV_DESEQ2_ANALYSIS${ params.assay_suffix }.csv" }

  label 'VV'

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)
    path(dge_table) // annotated dge table
    path(dge_table_rrnarm) // (rrna rm) annotated dge table

  output:
    path("VV_log.csv"), optional: params.skip_vv, emit: log

  script:
  def mode_params = params.mode == "microbes" ? "--mode microbes" : ""
  """
  if ${ !params.skip_vv } ; then
    vv_dge_deseq2.py --runsheet ${runsheet} --outdir ${publishdir} --assay_suffix ${params.assay_suffix} ${mode_params}
  fi
  """
}

process VV_STAR_ALIGNMENT {
  // Log publishing
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_VV_ALIGNMENT${ params.assay_suffix }.csv" }

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
    """
    # Run V&V unless user requests to skip V&V
    if ${ !params.skip_vv } ; then
      vv_star_alignment.py --runsheet ${runsheet} --outdir ${publishdir} --assay-suffix ${params.assay_suffix}
    fi
    """
}

process VV_RSEM_COUNTS {
  // Log publishing
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_VV_COUNTS${ params.assay_suffix }.csv" }

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
  
  script:
    """
    # Run V&V unless user requests to skip V&V
    if ${ !params.skip_vv } ; then
      vv_rsem_counts.py --runsheet ${runsheet} --outdir ${publishdir} --assay-suffix ${params.assay_suffix}
    fi
    """
}

process VV_CONCAT_FILTER {
  publishDir "${publishdir}/VV_Logs",
    mode: params.publish_dir_mode

  label 'VV'

  input:
    val(publishdir)
    path("VV_in.csv")

  output:
    tuple path("VV_log_final${params.assay_suffix}.csv"), path("VV_log_final_only_issues${params.assay_suffix}.csv")

  script:
    """
    concat_logs.py --assay_suffix ${params.assay_suffix}
    filter_to_only_issues.py --assay_suffix ${params.assay_suffix}
    """
}