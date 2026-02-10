include { RAW_READS_WORKFLOW; TRIMMED_READS_WORKFLOW; BAM_FILES_WORKFLOW; GENES_RESULTS_WORKFLOW; COUNTS_TABLE_WORKFLOW; DGE_TABLE_WORKFLOW } from '../subworkflows/rnaseq_subworkflows.nf'

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
                println "Executing TRIMMED_READS_WORKFLOW"
                TRIMMED_READS_WORKFLOW(
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
                
            case "bam_files":
                println "Executing BAM_FILES_WORKFLOW"
                BAM_FILES_WORKFLOW(
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
                
            case "genes_results":
                println "Executing GENES_RESULTS_WORKFLOW"
                GENES_RESULTS_WORKFLOW(
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
                
            case "counts_table":
                println "Executing COUNTS_TABLE_WORKFLOW"
                COUNTS_TABLE_WORKFLOW(
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

            case "dge_table":
                println "Executing DGE_TABLE_WORKFLOW"
                DGE_TABLE_WORKFLOW(
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
        }
}