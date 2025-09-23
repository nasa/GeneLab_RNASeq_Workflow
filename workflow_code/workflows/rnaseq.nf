include { RAW_READS_WORKFLOW } from '../subworkflows/rnaseq_subworkflows.nf'

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

workflow RNASEQ {
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
        println "Entry point: '${params.entry_point}'"
        
        switch(params.entry_point) {
            case "raw_reads":
                println "Executing RAW_READS_WORKFLOW"
                RAW_READS_WORKFLOW(
                    ch_outdir,
                    dp_tools_plugin,
                    annotations_csv_url_string,
                    accession,
                    isa_archive_path,
                    runsheet_path,
                    api_url,
                    reference_store_path,
                    derived_store_path,
                    reference_source,
                    reference_version,
                    reference_fasta,
                    reference_gtf
                )
                break
                
            case "trimmed_reads":
                println "TRIMMED_READS entry point not yet implemented"
                exit 1
                break
                
            case "bam_files":
                println "BAM_FILES entry point not yet implemented"
                exit 1
                break
                
            case "genes_results":
                println "GENES_RESULTS entry point not yet implemented"
                exit 1
                break
                
            case "counts_table":
                println "COUNTS_TABLE entry point not yet implemented"
                exit 1
                break
                
            case "dge_table":
                println "DGE_TABLE entry point not yet implemented"
                exit 1
                break
                
            default:
                println "ERROR: Invalid entry point '${params.entry_point}'"
                exit 1
        }
}