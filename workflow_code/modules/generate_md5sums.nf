process GENERATE_MD5SUMS {
  // Generates tabular data for GeneLab raw and processed data
    publishDir "${ch_outdir}/GeneLab",
        mode: params.publish_dir_mode,
        pattern: "*md5sum*"

    input:
        path(ch_outdir)
    output:
        path("raw_md5sum${params.assay_suffix}.tsv"), emit: raw_md5sum
        path("processed_md5sum${params.assay_suffix}.tsv"), emit: processed_md5sum

    script:
        def assay_suffix_arg = params.assay_suffix ? "--assay_suffix ${params.assay_suffix}" : ""
        """
        generate_md5sums.py --outdir ${ch_outdir} ${assay_suffix_arg}
        """
}