include { PARSE_RUNSHEET; mutate_to_single_end; PARSE_TRIMMED_RUNSHEET; PARSE_BAM_RUNSHEET } from './parse_runsheet.nf'
include { FETCH_ISA } from '../modules/fetch_isa.nf'
include { ISA_TO_RUNSHEET } from '../modules/isa_to_runsheet.nf'
include { GET_ACCESSIONS } from '../modules/get_accessions.nf'
include { STAGE_RAW_READS } from './stage_raw_reads.nf'
include { DOWNLOAD_OSDR_READS } from '../modules/download_osdr_reads.nf'
include { DOWNLOAD_OSDR_BAM } from '../modules/download_osdr_bam.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include { COPY_READS } from '../modules/copy_reads.nf'
include { COPY_BAMS } from '../modules/copy_bams.nf'


/**
 * For use within STAGE_ENTRY_TRIMMED_READS
 * 
 * 
 */
workflow STAGE_TRIMMED_READS {
    take:
        ch_outdir
        ch_samples

    main:
        truncate_to = params.truncate_to
        if ( truncate_to ) {
        // Download truncated raw reads
            ch_samples | map { it -> it[0].paired_end ? [it[0], it[1][0], it[1][1]] : [it[0], it[1][0]]}
                 | branch {
                   paired: it.size() == 3
                   single: it.size() == 2
                 }
                 | set{ ch_raw_read_pointers }
            
            
            // TO DO: Move the two splitFastq truncation steps into processes that run in parallel. Low priority since this is just for debugging.

            // PAIRED END
            // Only difference is the splitFastq arg 'pe'
            ch_raw_read_pointers.paired | splitFastq(pe: true, decompress: true, compress: true, limit: truncate_to, by: truncate_to, file: true)
                                        | map { it -> [ it[0], [ it[1], it[2] ] ]}
                                        // | view { it -> "TRUNCATED PAIRED READS ($truncate_to): $it[0]"}
                                        | set { ch_trimmed_reads }
            // SINGLE END
            // Only difference is the splitFastq arg 'pe'
            ch_raw_read_pointers.single | splitFastq(decompress: true, compress: true, limit: truncate_to, by: truncate_to, file: true)
                                        | map { it -> [ it[0], [ it[1] ] ]}
                                        // | view { it -> "TRUNCATED SINGLE READS ($truncate_to): $it[0]"}
                                        | mix( ch_trimmed_reads )
                                        | set { ch_trimmed_reads }

            // Rename and publish raw reads
            COPY_READS(ch_outdir, ch_trimmed_reads, "trimmed")
            // Collect sample IDs into a file
            COPY_READS.out.reads | map{ it -> it[1] } | collect | set { ch_all_trimmed_reads }
            COPY_READS.out.reads | map { it[0].id }
                            | collectFile(name: "samples.txt", sort: true, newLine: true)
                            | set { samples_txt }
        } else {
        // 
            ch_samples | map { it -> it[0].paired_end ? [it[0], [ it[1][0], it[1][1] ]] : [it[0], [it[1][0]]]}
                         | set { ch_trimmed_reads }

            /// Rename and publish trimmed reads
            COPY_READS(ch_outdir, ch_trimmed_reads, "trimmed")
            // Collect sample IDs into a file
            COPY_READS.out.reads | map{ it -> it[1] } | collect | set { ch_all_trimmed_reads }
            COPY_READS.out.reads | map { it[0].id }
                            | collectFile(name: "samples.txt", sort: true, newLine: true)
                            | set { samples_txt }
        }


    emit:
        trimmed_reads = COPY_READS.out.reads 
        ch_all_trimmed_reads = ch_all_trimmed_reads
        samples_txt =  samples_txt
}

/**
 * STAGE_ENTRY_TRIMMED_READS
 * 
 * This subworkflow handles the initial setup of the RNAseq analysis starting from trimmed reads:
 * 1. Sets up the output directory structure
 * 2. Fetches accessions if needed
 * 3. Obtains or creates the runsheet
 * 4. Parses the runsheet and stages trimmed reads
 */
workflow STAGE_ENTRY_TRIMMED_READS {
    take:
        ch_outdir
        dp_tools_plugin
        accession
        isa_archive_path
        runsheet_path
        api_url

    main:
        // Parse accession, structure output directory as:
        // params.outdir/
        //   ├── [GLDS-#|results]/ # Main pipeline results
        //   └── nextflow_info/    # Pipeline execution metadata
        Channel.empty() | set { osd_accession }
        Channel.empty() | set { glds_accession }
        
        if ( accession ) {
            GET_ACCESSIONS( accession, api_url )
            osd_accession = GET_ACCESSIONS.out.accessions_txt.map { it.readLines()[0].trim() }
            glds_accession = GET_ACCESSIONS.out.accessions_txt.map { it.readLines()[1].trim() }
            ch_outdir = ch_outdir.combine(glds_accession).map { outdir, glds -> "$outdir/$glds" }
        }
        else {
            ch_outdir = ch_outdir.map { it + "/results" }
        }
        ch_outdir = ch_outdir.first()

        Channel.empty() | set { isa_archive }
        if ( runsheet_path == null ) { // if runsheet_path is not provided, set it up from ISA input
            if ( isa_archive_path == null ) { // if isa_archive_path is not provided, fetch the ISA
                FETCH_ISA( ch_outdir, osd_accession, glds_accession )
                isa_archive = FETCH_ISA.out.isa_archive
            } else {
                // isa_archive_path is already a channel, use it directly
                isa_archive = isa_archive_path
            }
            ISA_TO_RUNSHEET( ch_outdir, osd_accession, glds_accession, isa_archive, dp_tools_plugin )
            runsheet_path = ISA_TO_RUNSHEET.out.runsheet

        } else if ( isa_archive_path != null ) {
            // if runsheet_path is provided and isa_archive_path is also provided, just pass through the provided ISA archive
            isa_archive = isa_archive_path
        }

        // Validate input parameters and runsheet
        validateParameters()
        
        // If entry point is trimmed_reads and no original runsheet was provided, download trimmed reads from OSDR 
        if ( params.entry_point == "trimmed_reads" && params.runsheet_path == null ) {
            PARSE_RUNSHEET( runsheet_path )
            samples = PARSE_RUNSHEET.out.samples
            DOWNLOAD_OSDR_READS( 
                ch_outdir, 
                osd_accession, 
                glds_accession,
                samples.map { meta, reads -> meta },
                "trimmed" 
            )
            
            // Check for download failures and exit gracefully if any occurred
            DOWNLOAD_OSDR_READS.out.failure_log
                | collect
                | subscribe { failure_files ->
                    if (failure_files.size() > 0) {
                        log.info "OSDR trimmed read file download failed for ${failure_files.size()} samples"
                        failure_files.each { file -> 
                            // Extract sample ID from filename like "FL_Bsu_1_GLbulkRNAseq_failed_download_trimmed.txt"
                            def sample_id = file.name.split("${params.assay_suffix}_failed_download_")[0]
                            log.info "OSDR trimmed read file download failed for ${sample_id}"
                        }
                        System.exit(0)  // Graceful exit without error
                    } else {
                        log.info "All trimmed read files downloaded successfully"
                    }
                }
            
            // Output is already in the right format: [meta, files]
            trimmed_reads = DOWNLOAD_OSDR_READS.out.trimmed_reads
 
        }

        else {
            PARSE_TRIMMED_RUNSHEET( runsheet_path )
            samples = PARSE_TRIMMED_RUNSHEET.out.samples
            STAGE_TRIMMED_READS( ch_outdir, samples )
            trimmed_reads = STAGE_TRIMMED_READS.out.trimmed_reads
        }
        samples_txt = trimmed_reads | map { it[0].id }
                            | collectFile(name: "samples.txt", sort: true, newLine: true) 


    emit:
        ch_outdir       = ch_outdir
        samples         = samples
        trimmed_reads   = trimmed_reads
        samples_txt     = samples_txt
        runsheet_path   = runsheet_path
        isa_archive     = isa_archive
        osd_accession   = osd_accession
        glds_accession  = glds_accession
}

/**
 * STAGE_ENTRY_BAM_FILES
 * 
 * This subworkflow handles the initial setup of the RNAseq analysis starting from BAM files:
 * 1. Sets up the output directory structure
 * 2. Fetches accessions if needed
 * 3. Obtains or creates the runsheet
 * 4. Stages BAM files (download from OSDR or use local paths)
 */
workflow STAGE_ENTRY_BAM_FILES {
    take:
        ch_outdir
        dp_tools_plugin
        accession
        isa_archive_path
        runsheet_path
        api_url

    main:
        // Parse accession, structure output directory as:
        // params.outdir/
        //   ├── [GLDS-#|results]/ # Main pipeline results
        //   └── nextflow_info/    # Pipeline execution metadata
        Channel.empty() | set { osd_accession }
        Channel.empty() | set { glds_accession }
        
        if ( accession ) {
            GET_ACCESSIONS( accession, api_url )
            osd_accession = GET_ACCESSIONS.out.accessions_txt.map { it.readLines()[0].trim() }
            glds_accession = GET_ACCESSIONS.out.accessions_txt.map { it.readLines()[1].trim() }
            ch_outdir = ch_outdir.combine(glds_accession).map { outdir, glds -> "$outdir/$glds" }
        }
        else {
            ch_outdir = ch_outdir.map { it + "/results" }
        }
        ch_outdir = ch_outdir.first()

        Channel.empty() | set { isa_archive }
        if ( runsheet_path == null ) {
            if ( isa_archive_path == null ) {
                FETCH_ISA( ch_outdir, osd_accession, glds_accession )
                isa_archive = FETCH_ISA.out.isa_archive
            } else {
                isa_archive = isa_archive_path
            }
            ISA_TO_RUNSHEET( ch_outdir, osd_accession, glds_accession, isa_archive, dp_tools_plugin )
            runsheet_path = ISA_TO_RUNSHEET.out.runsheet

        } else if ( isa_archive_path != null ) {
            isa_archive = isa_archive_path
        }

        // Validate input parameters and runsheet
        validateParameters()
        
        // If entry point is bam_files and no original runsheet was provided, download BAM files from OSDR 
        if ( params.entry_point == "bam_files" && params.runsheet_path == null ) {
            PARSE_RUNSHEET( runsheet_path )
            samples = PARSE_RUNSHEET.out.samples
            
            DOWNLOAD_OSDR_BAM( 
                ch_outdir, 
                osd_accession, 
                glds_accession,
                samples.map { meta, reads -> meta }
            )
            
            // Check for download failures and exit gracefully if any occurred
            DOWNLOAD_OSDR_BAM.out.failure_log
                | collect
                | subscribe { failure_files ->
                    if (failure_files.size() > 0) {
                        log.info "OSDR BAM file download failed for ${failure_files.size()} samples"
                        failure_files.each { file -> 
                            def sample_id = file.name.split("${params.assay_suffix}_failed_download_")[0]
                            log.info "OSDR BAM file download failed for ${sample_id}"
                        }
                        log.info "OSDR BAM file download failed for ${sample_id}"
                        System.exit(0)  // Graceful exit without error
                    } else {
                        log.info "All BAM files downloaded successfully"
                    }
                }
            
            bam_files = DOWNLOAD_OSDR_BAM.out.bam_files

        } else {
            PARSE_BAM_RUNSHEET( runsheet_path )
            samples = PARSE_BAM_RUNSHEET.out.samples
            runsheet_path = PARSE_BAM_RUNSHEET.out.runsheet
            
            // Rename BAM files to standard names
            COPY_BAMS(ch_outdir, samples)
            bam_files = COPY_BAMS.out.bam_files
            samples_txt = bam_files | map { it[0].id }
                            | collectFile(name: "samples.txt", sort: true, newLine: true)
        }

    emit:
        ch_outdir       = ch_outdir
        samples         = samples
        bam_files       = bam_files
        samples_txt     = samples_txt
        runsheet_path   = runsheet_path
        isa_archive     = isa_archive
        osd_accession   = osd_accession
        glds_accession  = glds_accession
} 