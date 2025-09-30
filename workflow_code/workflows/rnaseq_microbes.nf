include { RAW_READS_MICROBES_WORKFLOW } from '../subworkflows/rnaseq_microbes_subworkflows.nf'
include { TRIMMED_READS_MICROBES_WORKFLOW } from '../subworkflows/rnaseq_microbes_subworkflows.nf'

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

workflow RNASEQ_MICROBES {
    take:
        ch_outdir
        dp_tools_plugin
        annotations_csv_url_string
        accession
        isa_archive_path
        runsheet_path
        api_url
        force_single_end
        truncate_to
        reference_source
        reference_version
        reference_fasta
        reference_gtf
        reference_store_path
        derived_store_path
        entry_point
    main:
        // Simple entry point checking
        println "Entry point: '${params.entry_point}'"
        
        // Route to appropriate workflow based on entry point
        switch(params.entry_point) {
            case "raw_reads":
                println "Executing RAW_READS_MICROBES_WORKFLOW"
                RAW_READS_MICROBES_WORKFLOW(
                    ch_outdir,
                    dp_tools_plugin,
                    annotations_csv_url_string,
                    accession,
                    isa_archive_path,
                    runsheet_path,
                    api_url,
                    force_single_end,
                    truncate_to,
                    reference_source,
                    reference_version,
                    reference_fasta,
                    reference_gtf,
                    reference_store_path,
                    derived_store_path
                )
                break
                
            case "trimmed_reads":
                println "Executing TRIMMED_READS_MICROBES_WORKFLOW"
                TRIMMED_READS_MICROBES_WORKFLOW(
                    ch_outdir,
                    dp_tools_plugin,
                    annotations_csv_url_string,
                    accession,
                    isa_archive_path,
                    runsheet_path,
                    api_url,
                    force_single_end,
                    truncate_to,
                    reference_source,
                    reference_version,
                    reference_fasta,
                    reference_gtf,
                    reference_store_path,
                    derived_store_path
                )
                break
                
            case "bam_files":
                println "BAM_FILES entry point not yet implemented for microbes"
                exit 1
                break
                
            case "counts_table":
                println "COUNTS_TABLE entry point not yet implemented for microbes"
                exit 1
                break
                
            case "dge_table":
                println "DGE_TABLE entry point not yet implemented for microbes"
                exit 1
                break
                
            default:
                println "ERROR: Invalid entry point '${params.entry_point}' for microbes workflow"
                exit 1
        }
}