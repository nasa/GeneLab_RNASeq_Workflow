// RAW_READS_WORKFLOW - Complete workflow starting from raw reads
include { STAGE_ANALYSIS } from './stage_analysis.nf'
include { PARSE_ANNOTATIONS_TABLE } from '../modules/parse_annotations_table.nf'
include { DOWNLOAD_REFERENCES } from '../modules/download_references.nf'
include { SUBSAMPLE_GENOME } from '../modules/subsample_genome.nf'
include { DOWNLOAD_ERCC } from '../modules/download_ercc.nf'
include { CONCAT_ERCC } from '../modules/concat_ercc.nf'
include { GTF_TO_PRED } from '../modules/gtf_to_pred.nf'
include { PRED_TO_BED } from '../modules/pred_to_bed.nf'
include { STAGE_RAW_READS } from './stage_raw_reads.nf'
include { FASTQC as RAW_FASTQC } from '../modules/fastqc.nf'
include { GET_MAX_READ_LENGTH } from '../modules/get_max_read_length.nf'
include { TRIMGALORE } from '../modules/trimgalore.nf'
include { FASTQC as TRIMMED_FASTQC } from '../modules/fastqc.nf'
include { BUILD_STAR_INDEX } from '../modules/build_star_index.nf'
include { ALIGN_STAR } from '../modules/align_star.nf'
include { SORT_AND_INDEX_BAM } from '../modules/sort_and_index_bam.nf'
include { INFER_EXPERIMENT } from '../modules/rseqc.nf'
include { GENEBODY_COVERAGE } from '../modules/rseqc.nf'
include { INNER_DISTANCE } from '../modules/rseqc.nf'
include { READ_DISTRIBUTION } from '../modules/rseqc.nf'
include { ASSESS_STRANDEDNESS } from '../modules/assess_strandedness.nf'
include { BUILD_RSEM_INDEX } from '../modules/build_rsem_index.nf'
include { QUANTIFY_STAR_GENES } from '../modules/quantify_star_genes.nf'
include { COUNT_ALIGNED } from '../modules/count_aligned.nf' 
include { EXTRACT_RRNA } from '../modules/extract_rrna.nf'
include { REMOVE_RRNA } from '../modules/remove_rrna.nf'
include { QUANTIFY_RSEM_GENES } from '../modules/quantify_rsem_genes.nf'
include { DGE_DESEQ2 } from '../modules/dge_deseq2.nf'
include { DGE_DESEQ2 as DGE_DESEQ2_RRNA_RM } from '../modules/dge_deseq2.nf'   
include { 
    MULTIQC as RAW_READS_MULTIQC 
    MULTIQC as TRIMMED_READS_MULTIQC 
    MULTIQC as ALIGN_MULTIQC 
    MULTIQC as GENEBODY_COVERAGE_MULTIQC 
    MULTIQC as INFER_EXPERIMENT_MULTIQC 
    MULTIQC as INNER_DISTANCE_MULTIQC 
    MULTIQC as READ_DISTRIBUTION_MULTIQC 
    MULTIQC as COUNT_MULTIQC
} from '../modules/multiqc.nf'
include { PARSE_QC_METRICS } from '../modules/parse_qc_metrics.nf'
include { VV_RAW_READS;
    VV_TRIMMED_READS;
    VV_STAR_ALIGNMENT;
    VV_RSEQC;
    VV_RSEM_COUNTS;
    VV_DGE_DESEQ2;
    VV_CONCAT_FILTER } from '../modules/vv.nf'
include { SOFTWARE_VERSIONS } from '../modules/software_versions.nf'
include { GENERATE_PROTOCOL } from '../modules/generate_protocol.nf'
include { STAGE_ENTRY_TRIMMED_READS; STAGE_ENTRY_BAM_FILES; STAGE_ENTRY_GENES_RESULTS } from './stage_entry_points.nf'

workflow RAW_READS_WORKFLOW {
    take:
        ch_outdir
        dp_tools_plugin
        annotations_csv_url_string
        accession
        isa_archive_path
        runsheet_path
        api_url
        reference_store_path
        derived_store_path
        reference_source
        reference_version
        reference_fasta
        reference_gtf

    main:
        // Stage analysis setup (directory structure, inputs, and raw reads)
        STAGE_ANALYSIS(
            ch_outdir,
            dp_tools_plugin,
            accession,
            isa_archive_path,
            runsheet_path,
            api_url
        )
        ch_outdir = STAGE_ANALYSIS.out.ch_outdir
        samples = STAGE_ANALYSIS.out.samples
        raw_reads = STAGE_ANALYSIS.out.raw_reads
        samples_txt = STAGE_ANALYSIS.out.samples_txt
        runsheet_path = STAGE_ANALYSIS.out.runsheet_path
        isa_archive = STAGE_ANALYSIS.out.isa_archive
        osd_accession = STAGE_ANALYSIS.out.osd_accession
        glds_accession = STAGE_ANALYSIS.out.glds_accession
        
        // Get dataset-wide metadata
        samples | first 
                | map { meta, reads -> meta }
                | set { ch_meta }

        ch_meta | map { meta -> meta.organism_sci }
        | set { organism_sci }

        PARSE_ANNOTATIONS_TABLE( annotations_csv_url_string, organism_sci )

        // Use reference input and gene annotations file workflow params if provided
        if ( params.reference_fasta && params.reference_gtf ) {
            genome_references_pre_subsample = Channel.fromPath([params.reference_fasta, params.reference_gtf], checkIfExists: true ).toList()
            Channel.value( params.reference_source ) | set { reference_source }
            Channel.value( params.reference_version ) | set { reference_version }
            Channel.value( params.reference_fasta ) | set { reference_fasta_url }
            Channel.value( params.reference_gtf ) | set { reference_gtf_url }
            Channel.value( params.gene_annotations_file ) | set { gene_annotations_url }
        } else{
            // Use annotations table to get reference inputs, organism-specific gene annotations file
            reference_source = PARSE_ANNOTATIONS_TABLE.out.reference_source
            reference_version = PARSE_ANNOTATIONS_TABLE.out.reference_version
            reference_fasta_url = PARSE_ANNOTATIONS_TABLE.out.reference_fasta_url
            reference_gtf_url = PARSE_ANNOTATIONS_TABLE.out.reference_gtf_url
            gene_annotations_url = PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url
        }

        DOWNLOAD_REFERENCES( reference_store_path, organism_sci, reference_source, reference_version, reference_fasta_url, reference_gtf_url )
        genome_references_pre_subsample = DOWNLOAD_REFERENCES.out.reference_files

        // Genomic region subsampling step is used only for debugging / testing 
        if ( params.genome_subsample ) {
            SUBSAMPLE_GENOME( derived_store_path, organism_sci, genome_references_pre_subsample, reference_source, reference_version )
            SUBSAMPLE_GENOME.out.build | flatten | toList | set { genome_references_pre_ercc }
        } else {
            genome_references_pre_subsample | flatten | toList | set { genome_references_pre_ercc }
        }

        // Add ERCC Fasta and GTF to genome files
        DOWNLOAD_ERCC( ch_meta.map { it.has_ercc }, reference_store_path ).ifEmpty([file("ERCC92.fa"), file("ERCC92.gtf")]) | set { ch_maybe_ercc_refs }
        CONCAT_ERCC( reference_store_path, organism_sci, reference_source, reference_version, genome_references_pre_ercc, ch_maybe_ercc_refs, ch_meta.map { it.has_ercc } )
        .ifEmpty { genome_references_pre_ercc.value }  | set { genome_references }
        
        // Convert GTF file to RSeQC-compatible BED file
        GTF_TO_PRED(
            derived_store_path,
            organism_sci,
            reference_source,
            reference_version,
            genome_references | map { it[1] }
        )
        PRED_TO_BED( 
            derived_store_path,
            organism_sci,
            reference_source,
            reference_version,
            GTF_TO_PRED.out.genome_pred
        )
        genome_bed = PRED_TO_BED.out.genome_bed

        // Run FastQC on raw reads  
        RAW_FASTQC( ch_outdir.map { it + "/00-RawData/FastQC_Reports" }, raw_reads )
        RAW_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] } // Collect the raw read fastqc zip files
        | flatten
        | collect // Collect all zip files into a single list
        | set { raw_fastqc_zip } // Create a channel with all zip files
        
        // Get the max read length by parsing the raw read fastqc zip files
        GET_MAX_READ_LENGTH( raw_fastqc_zip )
        max_read_length = GET_MAX_READ_LENGTH.out.length | map { it.toString().toInteger() }
        
        // Trim raw reads
        TRIMGALORE( ch_outdir.map { it + "/01-TG_Preproc" }, raw_reads )
        trimmed_reads = TRIMGALORE.out.reads
        trimgalore_reports = TRIMGALORE.out.reports | collect

        // Run FastQC on trimmed reads
        TRIMMED_FASTQC( ch_outdir.map { it + "/01-TG_Preproc/FastQC_Reports" }, trimmed_reads )
        TRIMMED_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] }
        | flatten 
        | collect
        | set { trimmed_fastqc_zip }

        // Build STAR genome index
        BUILD_STAR_INDEX(derived_store_path, organism_sci, reference_source, reference_version, genome_references, ch_meta, max_read_length )
        star_index_dir = BUILD_STAR_INDEX.out.index_dir

        // STAR two-pass alignment
        ALIGN_STAR( ch_outdir.map { it + "/02-STAR_Alignment" }, trimmed_reads, star_index_dir )
        star_alignment_logs = ALIGN_STAR.out.alignment_logs | collect
        
        // Sort and index bam files
        SORT_AND_INDEX_BAM( ch_outdir.map { it + "/02-STAR_Alignment" }, ALIGN_STAR.out.bam_by_coord )
        sorted_bam = SORT_AND_INDEX_BAM.out.sorted_bam

        // RSeQC modules
        GENEBODY_COVERAGE( ch_outdir.map { it + "/RSeQC_Analyses/02_geneBody_coverage" }, sorted_bam, genome_bed )
        INFER_EXPERIMENT( ch_outdir.map { it + "/RSeQC_Analyses/03_infer_experiment" }, sorted_bam, genome_bed )
        INNER_DISTANCE( ch_outdir.map { it + "/RSeQC_Analyses/04_inner_distance" }, sorted_bam, genome_bed, max_read_length )
        READ_DISTRIBUTION( ch_outdir.map { it + "/RSeQC_Analyses/05_read_distribution" }, sorted_bam, genome_bed )
        infer_expt_out = INFER_EXPERIMENT.out.log | map { it[1] }
        | collect

        // Combine RSeQC module logs
        ch_rseqc_logs = Channel.empty()
        ch_rseqc_logs 
        | mix( INFER_EXPERIMENT.out.log_only,
                GENEBODY_COVERAGE.out.all_output,
                INNER_DISTANCE.out.all_output,
                READ_DISTRIBUTION.out.log_only )
                | collect
                | set{ ch_rseqc_logs }

        // Parse RSeQC infer_experiment.py results using thresholds set in bin/assess_strandedness.py to determine the strandedness of the dataset
        ASSESS_STRANDEDNESS( infer_expt_out )
        strandedness = ASSESS_STRANDEDNESS.out | map { it.text.split(":")[0] }

        // Create STAR counts table, nonzero gene counts
        QUANTIFY_STAR_GENES( ch_outdir.map { it + "/02-STAR_Alignment" }, samples_txt, ALIGN_STAR.out.reads_per_gene | toSortedList, strandedness )

        // Build RSEM transcriptome index
        BUILD_RSEM_INDEX(derived_store_path, organism_sci, reference_source, reference_version, genome_references, ch_meta )
        rsem_index_dir = BUILD_RSEM_INDEX.out.index_dir

        // Run RSEM on the transcriptome-aligned BAMs from STAR to calculate isoform-level transcript expression estimates and create a gene counts table
        COUNT_ALIGNED( ch_outdir.map { it + "/03-RSEM_Counts" }, ALIGN_STAR.out.bam_to_transcriptome, rsem_index_dir, strandedness )
        rsem_counts = COUNT_ALIGNED.out.counts | map { it[1] } | collect
        QUANTIFY_RSEM_GENES( ch_outdir.map { it + "/03-RSEM_Counts" }, samples_txt, rsem_counts )

        EXTRACT_RRNA ( organism_sci, genome_references | map { it[1] })
        REMOVE_RRNA ( ch_outdir.map { it + "/03-RSEM_Counts" }, EXTRACT_RRNA.out.rrna_ids, COUNT_ALIGNED.out.genes_results )

        dge_script = "${projectDir}/bin/dge_deseq2.Rmd"
        
        // Normalize counts, DGE, Add annotations to DGE table
        DGE_DESEQ2( ch_outdir, ch_meta, PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url, runsheet_path, COUNT_ALIGNED.out.genes_results.map{ it[1] } | collect, dge_script, "" )
        // For rRNArm counts: Normalize counts, DGE, Add annotations to DGE table
        DGE_DESEQ2_RRNA_RM( ch_outdir, ch_meta, PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url, runsheet_path, REMOVE_RRNA.out.genes_results_rrnarm | toSortedList, dge_script, "_rRNArm" )

        // MultiQC
        ch_multiqc_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config ) : Channel.fromPath("NO_FILE")
        RAW_READS_MULTIQC( ch_outdir.map { it + "/00-RawData/MultiQC_Reports" }, samples_txt, raw_fastqc_zip, ch_multiqc_config, "raw_")
        TRIMMED_READS_MULTIQC( ch_outdir.map { it + "/01-TG_Preproc/MultiQC_Reports" }, samples_txt, trimmed_fastqc_zip | concat( TRIMGALORE.out.reports ) | collect, ch_multiqc_config, "trimmed_")
        ALIGN_MULTIQC( ch_outdir.map { it + "/02-STAR_Alignment/MultiQC_Reports" }, samples_txt, star_alignment_logs, ch_multiqc_config, "align_")
        INFER_EXPERIMENT_MULTIQC( ch_outdir.map { it + "/RSeQC_Analyses/MultiQC_Reports" }, samples_txt, INFER_EXPERIMENT.out.log | map { it[1] } | collect, ch_multiqc_config, "infer_exp_")
        GENEBODY_COVERAGE_MULTIQC( ch_outdir.map { it + "/RSeQC_Analyses/MultiQC_Reports" }, samples_txt, GENEBODY_COVERAGE.out.log | map { it[1] } | collect, ch_multiqc_config, "geneBody_cov_")
        INNER_DISTANCE_MULTIQC( ch_outdir.map { it + "/RSeQC_Analyses/MultiQC_Reports" }, samples_txt, INNER_DISTANCE.out.log | map { it[1] } | collect, ch_multiqc_config, "inner_dist_")
        READ_DISTRIBUTION_MULTIQC( ch_outdir.map { it + "/RSeQC_Analyses/MultiQC_Reports" }, samples_txt, READ_DISTRIBUTION.out.log | map { it[1] } | collect, ch_multiqc_config, "read_dist_")
        COUNT_MULTIQC( ch_outdir.map { it + "/03-RSEM_Counts/MultiQC_Reports" }, samples_txt, rsem_counts, ch_multiqc_config, "RSEM_count_")

        // all_multiqc_input = raw_fastqc_zip
        //             | concat( trimgalore_reports )
        //             | concat( trimmed_fastqc_zip )
        //             | concat( star_alignment_logs )
        //             | concat( INFER_EXPERIMENT.out.log | map { it[1] } | collect )
        //             | concat( GENEBODY_COVERAGE.out.log | map { it[1] } | collect )
        //             | concat( INNER_DISTANCE.out.log | map { it[1] } | collect )
        //             | concat( READ_DISTRIBUTION.out.log | map { it[1] } | collect )
        //             | concat( rsem_counts )
        //             | collect
        // ALL_MULTIQC( ch_outdir.map { it + "/GeneLab" }, samples_txt, all_multiqc_input, ch_multiqc_config, "all_")

        // Parse QC metrics
        all_multiqc_output = RAW_READS_MULTIQC.out.data
            | concat( TRIMMED_READS_MULTIQC.out.data )
            | concat( ALIGN_MULTIQC.out.data )
            | concat( GENEBODY_COVERAGE_MULTIQC.out.data )
            | concat( INFER_EXPERIMENT_MULTIQC.out.data )
            | concat( INNER_DISTANCE_MULTIQC.out.data )
            | concat( READ_DISTRIBUTION_MULTIQC.out.data )
            | concat( COUNT_MULTIQC.out.data )
            | collect
        PARSE_QC_METRICS(
            ch_outdir,
            osd_accession,
            ch_meta,
            isa_archive.ifEmpty(file("ISA.zip")),  // Use a placeholder if isa_archive is empty
            all_multiqc_output,
            QUANTIFY_RSEM_GENES.out.publishables,
            runsheet_path
        )

        VV_RAW_READS(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            RAW_READS_MULTIQC.out.zipped_data
        )
        VV_TRIMMED_READS(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            TRIMMED_READS_MULTIQC.out.zipped_data
        )
        VV_STAR_ALIGNMENT(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            ALIGN_MULTIQC.out.zipped_data,
            QUANTIFY_STAR_GENES.out.publishables,
            SORT_AND_INDEX_BAM.out.bam_only_files | collect,
        )
        VV_RSEQC(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            GENEBODY_COVERAGE_MULTIQC.out.zipped_data,
            INFER_EXPERIMENT_MULTIQC.out.zipped_data,
            Channel.empty() | mix(INNER_DISTANCE_MULTIQC.out.zipped_data) | collect | ifEmpty({ file("PLACEHOLDER") }),
            READ_DISTRIBUTION_MULTIQC.out.zipped_data
        )
        VV_RSEM_COUNTS(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            COUNT_MULTIQC.out.zipped_data,
            REMOVE_RRNA.out.genes_results_rrnarm | collect,
            QUANTIFY_RSEM_GENES.out.publishables
        )
        VV_DGE_DESEQ2(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            DGE_DESEQ2.out.dge_table,
            DGE_DESEQ2_RRNA_RM.out.dge_table
        )
        // Concatenate and filter V&V logs
        VV_CONCAT_FILTER(
            ch_outdir,
            VV_RAW_READS.out.log
                | mix( 
                    VV_TRIMMED_READS.out.log,
                    VV_STAR_ALIGNMENT.out.log,
                    VV_RSEQC.out.log,
                    VV_RSEM_COUNTS.out.log,
                    VV_DGE_DESEQ2.out.log
                )
                | collect
        )

        // Software Version Capturing
        nf_version = '"NEXTFLOW":\n    nextflow: '.concat("${nextflow.version}\n")
        ch_nextflow_version = Channel.value(nf_version)
        ch_software_versions = Channel.empty()
        // Mix in versions from each process
        ch_software_versions = ch_software_versions
            | mix(GTF_TO_PRED.out.versions)
            | mix(PRED_TO_BED.out.versions)
            | mix(RAW_FASTQC.out.versions)
            | mix(TRIMGALORE.out.versions)
            | mix(ALIGN_STAR.out.versions)
            | mix(SORT_AND_INDEX_BAM.out.versions)
            | mix(INFER_EXPERIMENT.out.versions)
            | mix(GENEBODY_COVERAGE.out.versions)
            | mix(INNER_DISTANCE.out.versions)
            | mix(READ_DISTRIBUTION.out.versions)
            | mix(COUNT_ALIGNED.out.versions)
            | mix(RAW_READS_MULTIQC.out.versions)
            | mix(DGE_DESEQ2.out.versions)
            | mix(VV_RAW_READS.out.versions)
            | mix(ch_nextflow_version)
        // Process the versions:
        ch_software_versions 
            | unique  
            | collectFile(
                newLine: true
            )
            | set { ch_final_software_versions }
        // Convert software versions combined yaml to markdown table
        SOFTWARE_VERSIONS(ch_outdir, ch_final_software_versions)

        GENERATE_PROTOCOL(ch_outdir,
            ch_meta,
            strandedness,
            SOFTWARE_VERSIONS.out.software_versions_yaml,
            reference_source,
            reference_version,
            genome_references_pre_ercc,
            runsheet_path
        )
}

workflow TRIMMED_READS_WORKFLOW {
        take:
        ch_outdir
        dp_tools_plugin
        annotations_csv_url_string
        accession
        isa_archive_path
        runsheet_path
        api_url
        reference_store_path
        derived_store_path
        reference_source
        reference_version
        reference_fasta
        reference_gtf

    main:
        // Stage analysis setup (directory structure, inputs, and trimmed reads)
        STAGE_ENTRY_TRIMMED_READS(
            ch_outdir,
            dp_tools_plugin,
            accession,
            isa_archive_path,
            runsheet_path,
            api_url
        )
        ch_outdir = STAGE_ENTRY_TRIMMED_READS.out.ch_outdir
        samples = STAGE_ENTRY_TRIMMED_READS.out.samples
        trimmed_reads = STAGE_ENTRY_TRIMMED_READS.out.trimmed_reads
        samples_txt = STAGE_ENTRY_TRIMMED_READS.out.samples_txt
        runsheet_path = STAGE_ENTRY_TRIMMED_READS.out.runsheet_path
        isa_archive = STAGE_ENTRY_TRIMMED_READS.out.isa_archive
        osd_accession = STAGE_ENTRY_TRIMMED_READS.out.osd_accession
        glds_accession = STAGE_ENTRY_TRIMMED_READS.out.glds_accession

        // Get dataset-wide metadata
        samples | first 
                | map { meta, reads -> meta }
                | set { ch_meta }

        ch_meta | map { meta -> meta.organism_sci }
        | set { organism_sci }

        PARSE_ANNOTATIONS_TABLE( annotations_csv_url_string, organism_sci )

        // Use reference input and gene annotations file workflow params if provided
        if ( params.reference_fasta && params.reference_gtf ) {
            genome_references_pre_subsample = Channel.fromPath([params.reference_fasta, params.reference_gtf], checkIfExists: true ).toList()
            Channel.value( params.reference_source ) | set { reference_source }
            Channel.value( params.reference_version ) | set { reference_version }
            Channel.value( params.reference_fasta ) | set { reference_fasta_url }
            Channel.value( params.reference_gtf ) | set { reference_gtf_url }
            Channel.value( params.gene_annotations_file ) | set { gene_annotations_url }
        } else{
            // Use annotations table to get reference inputs, organism-specific gene annotations file
            reference_source = PARSE_ANNOTATIONS_TABLE.out.reference_source
            reference_version = PARSE_ANNOTATIONS_TABLE.out.reference_version
            reference_fasta_url = PARSE_ANNOTATIONS_TABLE.out.reference_fasta_url
            reference_gtf_url = PARSE_ANNOTATIONS_TABLE.out.reference_gtf_url
            gene_annotations_url = PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url
        }

        DOWNLOAD_REFERENCES( reference_store_path, organism_sci, reference_source, reference_version, reference_fasta_url, reference_gtf_url )
        genome_references_pre_subsample = DOWNLOAD_REFERENCES.out.reference_files

        // Genomic region subsampling step is used only for debugging / testing 
        if ( params.genome_subsample ) {
            SUBSAMPLE_GENOME( derived_store_path, organism_sci, genome_references_pre_subsample, reference_source, reference_version )
            SUBSAMPLE_GENOME.out.build | flatten | toList | set { genome_references_pre_ercc }
        } else {
            genome_references_pre_subsample | flatten | toList | set { genome_references_pre_ercc }
        }

        // Add ERCC Fasta and GTF to genome files
        DOWNLOAD_ERCC( ch_meta.map { it.has_ercc }, reference_store_path ).ifEmpty([file("ERCC92.fa"), file("ERCC92.gtf")]) | set { ch_maybe_ercc_refs }
        CONCAT_ERCC( reference_store_path, organism_sci, reference_source, reference_version, genome_references_pre_ercc, ch_maybe_ercc_refs, ch_meta.map { it.has_ercc } )
        .ifEmpty { genome_references_pre_ercc.value }  | set { genome_references }
        
        // Convert GTF file to RSeQC-compatible BED file
        GTF_TO_PRED(
            derived_store_path,
            organism_sci,
            reference_source,
            reference_version,
            genome_references | map { it[1] }
        )
        PRED_TO_BED( 
            derived_store_path,
            organism_sci,
            reference_source,
            reference_version,
            GTF_TO_PRED.out.genome_pred
        )
        genome_bed = PRED_TO_BED.out.genome_bed

        // Run FastQC on trimmed reads
        TRIMMED_FASTQC( ch_outdir.map { it + "/01-TG_Preproc/FastQC_Reports" }, trimmed_reads )
        TRIMMED_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] }
        | flatten 
        | collect
        | set { trimmed_fastqc_zip }

        // Get the max read length by parsing the trimmed read fastqc zip files
        GET_MAX_READ_LENGTH( trimmed_fastqc_zip )
        max_read_length = GET_MAX_READ_LENGTH.out.length | map { it.toString().toInteger() }

        // Build STAR genome index
        BUILD_STAR_INDEX(derived_store_path, organism_sci, reference_source, reference_version, genome_references, ch_meta, max_read_length )
        star_index_dir = BUILD_STAR_INDEX.out.index_dir

        // STAR two-pass alignment
        ALIGN_STAR( ch_outdir.map { it + "/02-STAR_Alignment" }, trimmed_reads, star_index_dir )
        star_alignment_logs = ALIGN_STAR.out.alignment_logs | collect
        
        // Sort and index bam files
        SORT_AND_INDEX_BAM( ch_outdir.map { it + "/02-STAR_Alignment" }, ALIGN_STAR.out.bam_by_coord )
        sorted_bam = SORT_AND_INDEX_BAM.out.sorted_bam

        // RSeQC modules
        GENEBODY_COVERAGE( ch_outdir.map { it + "/RSeQC_Analyses/02_geneBody_coverage" }, sorted_bam, genome_bed )
        INFER_EXPERIMENT( ch_outdir.map { it + "/RSeQC_Analyses/03_infer_experiment" }, sorted_bam, genome_bed )
        INNER_DISTANCE( ch_outdir.map { it + "/RSeQC_Analyses/04_inner_distance" }, sorted_bam, genome_bed, max_read_length )
        READ_DISTRIBUTION( ch_outdir.map { it + "/RSeQC_Analyses/05_read_distribution" }, sorted_bam, genome_bed )
        infer_expt_out = INFER_EXPERIMENT.out.log | map { it[1] }
        | collect

        // Combine RSeQC module logs
        ch_rseqc_logs = Channel.empty()
        ch_rseqc_logs 
        | mix( INFER_EXPERIMENT.out.log_only,
                GENEBODY_COVERAGE.out.all_output,
                INNER_DISTANCE.out.all_output,
                READ_DISTRIBUTION.out.log_only )
                | collect
                | set{ ch_rseqc_logs }

        // Parse RSeQC infer_experiment.py results using thresholds set in bin/assess_strandedness.py to determine the strandedness of the dataset
        ASSESS_STRANDEDNESS( infer_expt_out )
        strandedness = ASSESS_STRANDEDNESS.out | map { it.text.split(":")[0] }

        // Create STAR counts table, nonzero gene counts
        QUANTIFY_STAR_GENES( ch_outdir.map { it + "/02-STAR_Alignment" }, samples_txt, ALIGN_STAR.out.reads_per_gene | toSortedList, strandedness )

        // Build RSEM transcriptome index
        BUILD_RSEM_INDEX(derived_store_path, organism_sci, reference_source, reference_version, genome_references, ch_meta )
        rsem_index_dir = BUILD_RSEM_INDEX.out.index_dir

        // Run RSEM on the transcriptome-aligned BAMs from STAR to calculate isoform-level transcript expression estimates and create a gene counts table
        COUNT_ALIGNED( ch_outdir.map { it + "/03-RSEM_Counts" }, ALIGN_STAR.out.bam_to_transcriptome, rsem_index_dir, strandedness )
        rsem_counts = COUNT_ALIGNED.out.counts | map { it[1] } | collect
        QUANTIFY_RSEM_GENES( ch_outdir.map { it + "/03-RSEM_Counts" }, samples_txt, rsem_counts )

        EXTRACT_RRNA ( organism_sci, genome_references | map { it[1] })
        REMOVE_RRNA ( ch_outdir.map { it + "/03-RSEM_Counts" }, EXTRACT_RRNA.out.rrna_ids, COUNT_ALIGNED.out.genes_results )

        dge_script = "${projectDir}/bin/dge_deseq2.Rmd"
        
        // Normalize counts, DGE, Add annotations to DGE table
        DGE_DESEQ2( ch_outdir, ch_meta, PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url, runsheet_path, COUNT_ALIGNED.out.genes_results.map{ it[1] } | collect, dge_script, "" )
        // For rRNArm counts: Normalize counts, DGE, Add annotations to DGE table
        DGE_DESEQ2_RRNA_RM( ch_outdir, ch_meta, PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url, runsheet_path, REMOVE_RRNA.out.genes_results_rrnarm | toSortedList, dge_script, "_rRNArm" )

        // MultiQC
        ch_multiqc_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config ) : Channel.fromPath("NO_FILE")
        TRIMMED_READS_MULTIQC( ch_outdir.map { it + "/01-TG_Preproc/MultiQC_Reports" }, samples_txt, trimmed_fastqc_zip, ch_multiqc_config, "trimmed_")
        ALIGN_MULTIQC( ch_outdir.map { it + "/02-STAR_Alignment/MultiQC_Reports" }, samples_txt, star_alignment_logs, ch_multiqc_config, "align_")
        INFER_EXPERIMENT_MULTIQC( ch_outdir.map { it + "/RSeQC_Analyses/MultiQC_Reports" }, samples_txt, INFER_EXPERIMENT.out.log | map { it[1] } | collect, ch_multiqc_config, "infer_exp_")
        GENEBODY_COVERAGE_MULTIQC( ch_outdir.map { it + "/RSeQC_Analyses/MultiQC_Reports" }, samples_txt, GENEBODY_COVERAGE.out.log | map { it[1] } | collect, ch_multiqc_config, "geneBody_cov_")
        INNER_DISTANCE_MULTIQC( ch_outdir.map { it + "/RSeQC_Analyses/MultiQC_Reports" }, samples_txt, INNER_DISTANCE.out.log | map { it[1] } | collect, ch_multiqc_config, "inner_dist_")
        READ_DISTRIBUTION_MULTIQC( ch_outdir.map { it + "/RSeQC_Analyses/MultiQC_Reports" }, samples_txt, READ_DISTRIBUTION.out.log | map { it[1] } | collect, ch_multiqc_config, "read_dist_")
        COUNT_MULTIQC( ch_outdir.map { it + "/03-RSEM_Counts/MultiQC_Reports" }, samples_txt, rsem_counts, ch_multiqc_config, "RSEM_count_")

        // Parse QC metrics
        all_multiqc_output = TRIMMED_READS_MULTIQC.out.data
            | concat( ALIGN_MULTIQC.out.data )
            | concat( GENEBODY_COVERAGE_MULTIQC.out.data )
            | concat( INFER_EXPERIMENT_MULTIQC.out.data )
            | concat( INNER_DISTANCE_MULTIQC.out.data )
            | concat( READ_DISTRIBUTION_MULTIQC.out.data )
            | concat( COUNT_MULTIQC.out.data )
            | collect
        PARSE_QC_METRICS(
            ch_outdir,
            osd_accession,
            ch_meta,
            isa_archive.ifEmpty(file("ISA.zip")),  // Use a placeholder if isa_archive is empty
            all_multiqc_output,
            QUANTIFY_RSEM_GENES.out.publishables,
            runsheet_path
        )
        VV_TRIMMED_READS(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            TRIMMED_READS_MULTIQC.out.zipped_data
        )
        VV_STAR_ALIGNMENT(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            ALIGN_MULTIQC.out.zipped_data,
            QUANTIFY_STAR_GENES.out.publishables,
            SORT_AND_INDEX_BAM.out.bam_only_files | collect,
        )
        VV_RSEQC(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            GENEBODY_COVERAGE_MULTIQC.out.zipped_data,
            INFER_EXPERIMENT_MULTIQC.out.zipped_data,
            Channel.empty() | mix(INNER_DISTANCE_MULTIQC.out.zipped_data) | collect | ifEmpty({ file("PLACEHOLDER") }),
            READ_DISTRIBUTION_MULTIQC.out.zipped_data
        )
        VV_RSEM_COUNTS(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            COUNT_MULTIQC.out.zipped_data,
            REMOVE_RRNA.out.genes_results_rrnarm | collect,
            QUANTIFY_RSEM_GENES.out.publishables
        )
        VV_DGE_DESEQ2(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            DGE_DESEQ2.out.dge_table,
            DGE_DESEQ2_RRNA_RM.out.dge_table
        )
        // Concatenate and filter V&V logs
        VV_CONCAT_FILTER(
            ch_outdir,
            VV_TRIMMED_READS.out.log
                | mix( 
                    VV_STAR_ALIGNMENT.out.log,
                    VV_RSEQC.out.log,
                    VV_RSEM_COUNTS.out.log,
                    VV_DGE_DESEQ2.out.log
                )
                | collect
        )

        // Software Version Capturing
        nf_version = '"NEXTFLOW":\n    nextflow: '.concat("${nextflow.version}\n")
        ch_nextflow_version = Channel.value(nf_version)
        ch_software_versions = Channel.empty()
        // Mix in versions from each process
        ch_software_versions = ch_software_versions
            | mix(GTF_TO_PRED.out.versions)
            | mix(PRED_TO_BED.out.versions)
            | mix(ALIGN_STAR.out.versions)
            | mix(SORT_AND_INDEX_BAM.out.versions)
            | mix(INFER_EXPERIMENT.out.versions)
            | mix(GENEBODY_COVERAGE.out.versions)
            | mix(INNER_DISTANCE.out.versions)
            | mix(READ_DISTRIBUTION.out.versions)
            | mix(COUNT_ALIGNED.out.versions)
            | mix(TRIMMED_READS_MULTIQC.out.versions)
            | mix(DGE_DESEQ2.out.versions)
            | mix(VV_TRIMMED_READS.out.versions)
            | mix(ch_nextflow_version)
        // Process the versions:
        ch_software_versions 
            | unique  
            | collectFile(
                newLine: true
            )
            | set { ch_final_software_versions }
        // Convert software versions combined yaml to markdown table
        SOFTWARE_VERSIONS(ch_outdir, ch_final_software_versions)

        GENERATE_PROTOCOL(ch_outdir,
            ch_meta,
            strandedness,
            SOFTWARE_VERSIONS.out.software_versions_yaml,
            reference_source,
            reference_version,
            genome_references_pre_ercc,
            runsheet_path
        )
}

workflow BAM_FILES_WORKFLOW {
    take:
        ch_outdir
        dp_tools_plugin
        annotations_csv_url_string
        accession
        isa_archive_path
        runsheet_path
        api_url
        reference_store_path
        derived_store_path
        reference_source
        reference_version
        reference_fasta
        reference_gtf

    main:
        // Stage analysis setup (directory structure, inputs, and bam files)
        STAGE_ENTRY_BAM_FILES(
            ch_outdir,
            dp_tools_plugin,
            accession,
            isa_archive_path,
            runsheet_path,
            api_url
        )
        ch_outdir = STAGE_ENTRY_BAM_FILES.out.ch_outdir
        samples = STAGE_ENTRY_BAM_FILES.out.samples
        bam_files = STAGE_ENTRY_BAM_FILES.out.bam_files
        samples_txt = STAGE_ENTRY_BAM_FILES.out.samples_txt
        runsheet_path = STAGE_ENTRY_BAM_FILES.out.runsheet_path
        isa_archive = STAGE_ENTRY_BAM_FILES.out.isa_archive
        osd_accession = STAGE_ENTRY_BAM_FILES.out.osd_accession
        glds_accession = STAGE_ENTRY_BAM_FILES.out.glds_accession

        // Get dataset-wide metadata
        samples | first | map { meta, reads -> meta } | set { ch_meta }
        ch_meta | map { meta -> meta.organism_sci } | set { organism_sci }

        PARSE_ANNOTATIONS_TABLE( annotations_csv_url_string, organism_sci )

        // Use reference input and gene annotations file workflow params if provided
        if ( params.reference_fasta && params.reference_gtf ) {
            genome_references_pre_subsample = Channel.fromPath([params.reference_fasta, params.reference_gtf], checkIfExists: true ).toList()
            Channel.value( params.reference_source ) | set { reference_source }
            Channel.value( params.reference_version ) | set { reference_version }
            Channel.value( params.reference_fasta ) | set { reference_fasta_url }
            Channel.value( params.reference_gtf ) | set { reference_gtf_url }
            Channel.value( params.gene_annotations_file ) | set { gene_annotations_url }
        } else{
            // Use annotations table to get reference inputs, organism-specific gene annotations file
            reference_source = PARSE_ANNOTATIONS_TABLE.out.reference_source
            reference_version = PARSE_ANNOTATIONS_TABLE.out.reference_version
            reference_fasta_url = PARSE_ANNOTATIONS_TABLE.out.reference_fasta_url
            reference_gtf_url = PARSE_ANNOTATIONS_TABLE.out.reference_gtf_url
            gene_annotations_url = PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url
        }

        DOWNLOAD_REFERENCES( reference_store_path, organism_sci, reference_source, reference_version, reference_fasta_url, reference_gtf_url )
        genome_references_pre_subsample = DOWNLOAD_REFERENCES.out.reference_files

        // Genomic region subsampling step is used only for debugging / testing 
        if ( params.genome_subsample ) {
            SUBSAMPLE_GENOME( derived_store_path, organism_sci, genome_references_pre_subsample, reference_source, reference_version )
            SUBSAMPLE_GENOME.out.build | flatten | toList | set { genome_references_pre_ercc }
        } else {
            genome_references_pre_subsample | flatten | toList | set { genome_references_pre_ercc }
        }


        // Add ERCC Fasta and GTF to genome files
        DOWNLOAD_ERCC( ch_meta.map { it.has_ercc }, reference_store_path ).ifEmpty([file("ERCC92.fa"), file("ERCC92.gtf")]) | set { ch_maybe_ercc_refs }
        CONCAT_ERCC( reference_store_path, organism_sci, reference_source, reference_version, genome_references_pre_ercc, ch_maybe_ercc_refs, ch_meta.map { it.has_ercc } )
        .ifEmpty { genome_references_pre_ercc.value }  | set { genome_references }
        
        // Convert GTF file to RSeQC-compatible BED file
        GTF_TO_PRED(
            derived_store_path,
            organism_sci,
            reference_source,
            reference_version,
            genome_references | map { it[1] }
        )
        PRED_TO_BED( 
            derived_store_path,
            organism_sci,
            reference_source,
            reference_version,
            GTF_TO_PRED.out.genome_pred
        )
        genome_bed = PRED_TO_BED.out.genome_bed
        
        // For BAM entry point: Convert params.strandedness to expected values
        def params_strandedness = params.strandedness ?: "none"
        def converted_strandedness = params_strandedness == "forward" ? "sense" : params_strandedness == "reverse" ? "antisense" : params_strandedness == "none" ? "unstranded" : params_strandedness
        strandedness = Channel.value(converted_strandedness)

        // Build RSEM transcriptome index
        BUILD_RSEM_INDEX(derived_store_path, organism_sci, reference_source, reference_version, genome_references, ch_meta )
        rsem_index_dir = BUILD_RSEM_INDEX.out.index_dir

        // Run RSEM on the transcriptome-aligned BAMs from STAR to calculate isoform-level transcript expression estimates and create a gene counts table
        COUNT_ALIGNED( ch_outdir.map { it + "/03-RSEM_Counts" }, bam_files, rsem_index_dir, strandedness )
        rsem_counts = COUNT_ALIGNED.out.counts | map { it[1] } | collect
        QUANTIFY_RSEM_GENES( ch_outdir.map { it + "/03-RSEM_Counts" }, samples_txt, rsem_counts )

        EXTRACT_RRNA ( organism_sci, genome_references | map { it[1] })
        REMOVE_RRNA ( ch_outdir.map { it + "/03-RSEM_Counts" }, EXTRACT_RRNA.out.rrna_ids, COUNT_ALIGNED.out.genes_results )

        dge_script = "${projectDir}/bin/dge_deseq2.Rmd"
        
        // Normalize counts, DGE, Add annotations to DGE table
        DGE_DESEQ2( ch_outdir, ch_meta, PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url, runsheet_path, COUNT_ALIGNED.out.genes_results.map{ it[1] } | collect, dge_script, "" )
        // For rRNArm counts: Normalize counts, DGE, Add annotations to DGE table
        DGE_DESEQ2_RRNA_RM( ch_outdir, ch_meta, PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url, runsheet_path, REMOVE_RRNA.out.genes_results_rrnarm | toSortedList, dge_script, "_rRNArm" )

        // MultiQC
        ch_multiqc_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config ) : Channel.fromPath("NO_FILE")
        COUNT_MULTIQC( ch_outdir.map { it + "/03-RSEM_Counts/MultiQC_Reports" }, samples_txt, rsem_counts, ch_multiqc_config, "RSEM_count_")

        // Parse QC metrics
        all_multiqc_output = COUNT_MULTIQC.out.data
        PARSE_QC_METRICS(
            ch_outdir,
            osd_accession,
            ch_meta,
            isa_archive.ifEmpty(file("ISA.zip")),  // Use a placeholder if isa_archive is empty
            all_multiqc_output,
            QUANTIFY_RSEM_GENES.out.publishables,
            runsheet_path
        )
        VV_RSEM_COUNTS(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            COUNT_MULTIQC.out.zipped_data,
            REMOVE_RRNA.out.genes_results_rrnarm | collect,
            QUANTIFY_RSEM_GENES.out.publishables
        )
        VV_DGE_DESEQ2(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            DGE_DESEQ2.out.dge_table,
            DGE_DESEQ2_RRNA_RM.out.dge_table
        )
        // Concatenate and filter V&V logs
        VV_CONCAT_FILTER(
            ch_outdir,
            VV_RSEM_COUNTS.out.log
                | mix( 
                    VV_DGE_DESEQ2.out.log
                )
                | collect
        )

        // Software Version Capturing
        nf_version = '"NEXTFLOW":\n    nextflow: '.concat("${nextflow.version}\n")
        ch_nextflow_version = Channel.value(nf_version)
        ch_software_versions = Channel.empty()
        // Mix in versions from each process
        ch_software_versions = ch_software_versions
            | mix(GTF_TO_PRED.out.versions)
            | mix(PRED_TO_BED.out.versions)
            | mix(COUNT_ALIGNED.out.versions)
            | mix(COUNT_MULTIQC.out.versions)
            | mix(DGE_DESEQ2.out.versions)
            | mix(VV_RSEM_COUNTS.out.versions)
            | mix(ch_nextflow_version)
        // Process the versions:
        ch_software_versions 
            | unique  
            | collectFile(
                newLine: true
            )
            | set { ch_final_software_versions }
        // Convert software versions combined yaml to markdown table
        SOFTWARE_VERSIONS(ch_outdir, ch_final_software_versions)

        GENERATE_PROTOCOL(ch_outdir,
            ch_meta,
            strandedness,
            SOFTWARE_VERSIONS.out.software_versions_yaml,
            reference_source,
            reference_version,
            genome_references_pre_ercc,
            runsheet_path
        )
}

workflow GENES_RESULTS_WORKFLOW {
    take:
        ch_outdir
        dp_tools_plugin
        annotations_csv_url_string
        accession
        isa_archive_path
        runsheet_path
        api_url
        reference_store_path
        derived_store_path
        reference_source
        reference_version
        reference_fasta
        reference_gtf

    main:
        // Stage analysis setup (directory structure, inputs, and RSEM genes.results files)
        STAGE_ENTRY_GENES_RESULTS(
            ch_outdir,
            dp_tools_plugin,
            accession,
            isa_archive_path,
            runsheet_path,
            api_url
        )
        ch_outdir = STAGE_ENTRY_GENES_RESULTS.out.ch_outdir
        samples = STAGE_ENTRY_GENES_RESULTS.out.samples
        genes_results = STAGE_ENTRY_GENES_RESULTS.out.genes_results
        samples_txt = STAGE_ENTRY_GENES_RESULTS.out.samples_txt
        runsheet_path = STAGE_ENTRY_GENES_RESULTS.out.runsheet_path
        isa_archive = STAGE_ENTRY_GENES_RESULTS.out.isa_archive
        osd_accession = STAGE_ENTRY_GENES_RESULTS.out.osd_accession
        glds_accession = STAGE_ENTRY_GENES_RESULTS.out.glds_accession

        // Get dataset-wide metadata
        samples | first | map { meta, reads -> meta } | set { ch_meta }
        ch_meta | map { meta -> meta.organism_sci } | set { organism_sci }

        PARSE_ANNOTATIONS_TABLE( annotations_csv_url_string, organism_sci )

        // Use reference input and gene annotations file workflow params if provided
        if ( params.reference_fasta && params.reference_gtf ) {
            genome_references_pre_subsample = Channel.fromPath([params.reference_fasta, params.reference_gtf], checkIfExists: true ).toList()
            Channel.value( params.reference_source ) | set { reference_source }
            Channel.value( params.reference_version ) | set { reference_version }
            Channel.value( params.reference_fasta ) | set { reference_fasta_url }
            Channel.value( params.reference_gtf ) | set { reference_gtf_url }
            Channel.value( params.gene_annotations_file ) | set { gene_annotations_url }
        } else{
            // Use annotations table to get reference inputs, organism-specific gene annotations file
            reference_source = PARSE_ANNOTATIONS_TABLE.out.reference_source
            reference_version = PARSE_ANNOTATIONS_TABLE.out.reference_version
            reference_fasta_url = PARSE_ANNOTATIONS_TABLE.out.reference_fasta_url
            reference_gtf_url = PARSE_ANNOTATIONS_TABLE.out.reference_gtf_url
            gene_annotations_url = PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url
        }

        DOWNLOAD_REFERENCES( reference_store_path, organism_sci, reference_source, reference_version, reference_fasta_url, reference_gtf_url )
        genome_references_pre_subsample = DOWNLOAD_REFERENCES.out.reference_files

        // Genomic region subsampling step is used only for debugging / testing 
        if ( params.genome_subsample ) {
            SUBSAMPLE_GENOME( derived_store_path, organism_sci, genome_references_pre_subsample, reference_source, reference_version )
            SUBSAMPLE_GENOME.out.build | flatten | toList | set { genome_references_pre_ercc }
        } else {
            genome_references_pre_subsample | flatten | toList | set { genome_references_pre_ercc }
        }


        // Add ERCC Fasta and GTF to genome files
        DOWNLOAD_ERCC( ch_meta.map { it.has_ercc }, reference_store_path ).ifEmpty([file("ERCC92.fa"), file("ERCC92.gtf")]) | set { ch_maybe_ercc_refs }
        CONCAT_ERCC( reference_store_path, organism_sci, reference_source, reference_version, genome_references_pre_ercc, ch_maybe_ercc_refs, ch_meta.map { it.has_ercc } )
        .ifEmpty { genome_references_pre_ercc.value }  | set { genome_references }
        
        // For genes results entry point: Convert params.strandedness to expected values
        def params_strandedness = params.strandedness ?: "none"
        def converted_strandedness = params_strandedness == "forward" ? "sense" : params_strandedness == "reverse" ? "antisense" : params_strandedness == "none" ? "unstranded" : params_strandedness
        strandedness = Channel.value(converted_strandedness)

        // Pass in genes results entry point files
        rsem_counts = genes_results | map { it[1] } | collect
        QUANTIFY_RSEM_GENES( ch_outdir.map { it + "/03-RSEM_Counts" }, samples_txt, rsem_counts )

        EXTRACT_RRNA ( organism_sci, genome_references | map { it[1] })
        REMOVE_RRNA ( ch_outdir.map { it + "/03-RSEM_Counts" }, EXTRACT_RRNA.out.rrna_ids, genes_results )

        dge_script = "${projectDir}/bin/dge_deseq2.Rmd"
        
        // Normalize counts, DGE, Add annotations to DGE table
        DGE_DESEQ2( ch_outdir, ch_meta, PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url, runsheet_path, rsem_counts, dge_script, "" )
        // For rRNArm counts: Normalize counts, DGE, Add annotations to DGE table
        DGE_DESEQ2_RRNA_RM( ch_outdir, ch_meta, PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url, runsheet_path, REMOVE_RRNA.out.genes_results_rrnarm | toSortedList, dge_script, "_rRNArm" )

        // MultiQC
        ch_multiqc_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config ) : Channel.fromPath("NO_FILE")

        // Parse QC metrics
        all_multiqc_output = Channel.empty()
        PARSE_QC_METRICS(
            ch_outdir,
            osd_accession,
            ch_meta,
            isa_archive.ifEmpty(file("ISA.zip")),  // Use a placeholder if isa_archive is empty
            all_multiqc_output,
            QUANTIFY_RSEM_GENES.out.publishables,
            runsheet_path
        )
        VV_DGE_DESEQ2(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            DGE_DESEQ2.out.dge_table,
            DGE_DESEQ2_RRNA_RM.out.dge_table
        )
        // Concatenate and filter V&V logs
        VV_CONCAT_FILTER(
            ch_outdir,
            VV_DGE_DESEQ2.out.log
        )

        // Software Version Capturing
        nf_version = '"NEXTFLOW":\n    nextflow: '.concat("${nextflow.version}\n")
        ch_nextflow_version = Channel.value(nf_version)
        ch_software_versions = Channel.empty()
        // Mix in versions from each process
        ch_software_versions = ch_software_versions
            | mix(DGE_DESEQ2.out.versions)
            | mix(VV_DGE_DESEQ2.out.versions)
            | mix(ch_nextflow_version)
        // Process the versions:
        ch_software_versions 
            | unique  
            | collectFile(
                newLine: true
            )
            | set { ch_final_software_versions }
        // Convert software versions combined yaml to markdown table
        SOFTWARE_VERSIONS(ch_outdir, ch_final_software_versions)

        GENERATE_PROTOCOL(ch_outdir,
            ch_meta,
            strandedness,
            SOFTWARE_VERSIONS.out.software_versions_yaml,
            reference_source,
            reference_version,
            genome_references_pre_ercc,
            runsheet_path
        )
}