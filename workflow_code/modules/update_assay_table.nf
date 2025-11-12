process UPDATE_ASSAY_TABLE {
    publishDir "${ch_processed_directory}/GeneLab/updated_curation_tables",
    mode: params.publish_dir_mode,
    pattern: "a_*.txt"

    input:
        path(ch_processed_directory)
    output:
        path("a_*.txt"), emit: assay_table, optional: true

    script:
    def mode_param = params.mode == "microbes" ? "--mode microbes" : ""
    def assay_suffix_arg = params.assay_suffix ? "--assay_suffix ${params.assay_suffix}" : ""
    """
    update_assay_table.py --outdir ${ch_processed_directory} ${mode_param} ${assay_suffix_arg} --glds_accession ${params.accession}
    """
}