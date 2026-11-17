include { STAGE } from '../subworkflows/stage.nf'

include { PARSE_ANNOTATIONS_TABLE } from '../modules/parse_annotations_table.nf'
include { FETCH_TABLE as FETCH_ANNOTATIONS_TABLE } from '../modules/fetch_remote.nf'
include { FETCH_TABLE as FETCH_ANNOTATIONS } from '../modules/fetch_remote.nf'
include { DOWNLOAD_REFERENCES; COPY_REFERENCES } from '../modules/download_references.nf'
include { SUBSAMPLE_GENOME } from '../modules/subsample_genome.nf'
include { DOWNLOAD_ERCC } from '../modules/download_ercc.nf'
include { CONCAT_ERCC } from '../modules/concat_ercc.nf'
include { ERCC_ANALYSIS } from '../modules/ercc_analysis.nf'
include { GTF_TO_PRED } from '../modules/gtf_to_pred.nf'
include { PRED_TO_BED } from '../modules/pred_to_bed.nf'

include { FASTQC as RAW_FASTQC } from '../modules/fastqc.nf'
include { FASTQC as TRIMMED_FASTQC } from '../modules/fastqc.nf'
include { GET_MAX_READ_LENGTH } from '../modules/get_max_read_length.nf'
include { TRIMGALORE } from '../modules/trimgalore.nf'

include { BUILD_STAR_INDEX } from '../modules/build_star_index.nf'
include { ALIGN_STAR } from '../modules/align_star.nf'
include { BUILD_RSEM_INDEX } from '../modules/build_rsem_index.nf'
include { QUANTIFY_STAR_GENES } from '../modules/quantify_star_genes.nf'
include { COUNT_ALIGNED } from '../modules/count_aligned.nf'
include { QUANTIFY_RSEM_GENES } from '../modules/quantify_rsem_genes.nf'

include { BUILD_BOWTIE2_INDEX } from '../modules/build_bowtie2_index.nf'
include { ALIGN_BOWTIE2 } from '../modules/align_bowtie2.nf'
include { SAM_TO_BAM } from '../modules/sam_to_bam.nf'
include { GET_GTF_FEATURES } from '../modules/get_gtf_features.nf'
include { FEATURECOUNTS } from '../modules/featurecounts.nf'
include { QUANTIFY_FEATURECOUNTS_GENES } from '../modules/quantify_featurecounts_genes.nf'

include { SORT_AND_INDEX_BAM } from '../modules/sort_and_index_bam.nf'

include { INFER_EXPERIMENT; GENEBODY_COVERAGE; INNER_DISTANCE; READ_DISTRIBUTION } from '../modules/rseqc.nf'
include { BAM_STRANDEDNESS } from '../subworkflows/bam_strandedness.nf'
include { ASSESS_STRANDEDNESS } from '../modules/assess_strandedness.nf'

include { EXTRACT_RRNA } from '../modules/extract_rrna.nf'
include { REMOVE_RRNA } from '../modules/remove_rrna.nf'
include { REMOVE_RRNA_FEATURECOUNTS } from '../modules/remove_rrna_featurecounts.nf'
include { REMOVE_RRNA_COUNTS_TABLE } from '../modules/remove_rrna_counts_table.nf'
include { REMOVE_RRNA_COUNTS_TABLE as REMOVE_RRNA_RSEM_TABLE } from '../modules/remove_rrna_counts_table.nf'
include { REMOVE_RRNA_COUNTS_TABLE as REMOVE_RRNA_FC_TABLE } from '../modules/remove_rrna_counts_table.nf'
include { DGE_DESEQ2 } from '../modules/dge_deseq2.nf'
include { DGE_DESEQ2 as DGE_DESEQ2_RRNA_RM } from '../modules/dge_deseq2.nf'
include { SPLIT_ORGANISM_PART_RUNSHEETS } from '../modules/split_organism_part_runsheets.nf'
include { FILTER_UNALIGNABLE } from '../modules/filter_unalignable.nf'
include { ANNOTATE_DGE_TABLE } from '../modules/annotate_dge_table.nf'

include { MULTIQC as RAW_READS_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as TRIMMED_READS_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as ALIGN_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as COUNT_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as GENEBODY_COVERAGE_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as INFER_EXPERIMENT_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as INNER_DISTANCE_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as READ_DISTRIBUTION_MULTIQC } from '../modules/multiqc.nf'

include { PARSE_QC_METRICS } from '../modules/parse_qc_metrics.nf'
include { VV_FASTQ as VV_RAW_FASTQ;
    VV_FASTQ as VV_TRIMMED_FASTQ;
    VV_RAW_READS;
    VV_TRIMMED_READS;
    VV_STAR_ALIGNMENT;
    VV_BOWTIE2_ALIGNMENT;
    VV_RSEQC;
    VV_RSEM_COUNTS;
    VV_FEATURECOUNTS;
    VV_DGE_DESEQ2;
    VV_CONCAT_FILTER } from '../modules/vv.nf'
include { SOFTWARE_VERSIONS } from '../modules/software_versions.nf'
include { GENERATE_PROTOCOL } from '../modules/generate_protocol.nf'

def strandedness_set() {
    return params.strandedness in ['none', 'forward', 'reverse']
}

def convert_strandedness(p) {
    if ( p == 'forward' ) { return 'sense' }
    if ( p == 'reverse' ) { return 'antisense' }
    return 'unstranded'
}

def as_list(x) {
    if ( x == null ) return []
    if ( x instanceof Collection && !(x instanceof CharSequence) ) return x as List
    return [x]
}

def dest_file(root, dir, f) {
    return ["${root}/${dir}/${f.name}".toString(), f]
}

def vv_log(root, name, f) {
    return ["${root}/VV_Logs/VV_log_${name}${params.assay_suffix}.csv".toString(), f]
}

def files_only(items) {
    return items.collectMany { x -> x instanceof Map ? [] : as_list(x) }
}

def flatten_fastqs(ch) {
    return ch.flatMap { row -> files_only(as_list(row)) }
}

def vv_fastq_reports(ch) {
    if ( params.skip_vv ) {
        return channel.fromPath("${projectDir}/assets/empty").collect()
    }
    return ch.collect()
}

def pub(ch, root, dir) {
    return ch.combine(root).combine(channel.value(dir)).flatMap { row ->
        def items = as_list(row)
        def d = items[-1]
        def r = items[-2]
        files_only(items[0..-3]).collect { f -> dest_file(r, d, f) }
    }
}

def pub_mqc(mqc, root, dir) {
    return pub(mqc.html, root, dir).mix(pub(mqc.zipped_data, root, dir))
}

def pub_by_meta(ch, root, dir) {
    return ch.combine(root).combine(channel.value(dir)).flatMap { row ->
        def items = as_list(row)
        def d = items[-1]
        def r = items[-2]
        def meta = items[0]
        files_only(items[1..-3]).collect { f -> dest_file(r, "${d}/${meta.id}", f) }
    }
}

def dge_jobs_from_manifest(manifest, runsheets) {
    return manifest.splitCsv(header: true, sep: '\t')
        .map { row -> [row.filename, row.output_label ?: ""] }
        .join(runsheets.flatten().map { f -> [f.name, f] })
        .map { row ->
            def items = as_list(row)
            [items[1], items[2]]
        }
}

def dge_input(jobs, meta, annot, counts, rrna) {
    return jobs.combine(meta).combine(annot).combine(counts).map { row ->
        def items = as_list(row)
        def lab = (items[0] ?: "").toString()
        def rs = items[1]
        def m = items[2]
        def a = items[3]
        def c = items.size() > 5 ? items[4..-1] : items[4]
        [m, a, rrna ? "${lab}_rRNArm".toString() : lab, rs, c]
    }
}

def pub_dge(dge, root, side_tables=true) {
    def normalized = dge.norm_counts.map { norm, _unnorm -> norm }
    def published = pub(normalized, root, '04-DESeq2_NormCounts')
        .mix(pub(dge.vst_norm_counts, root, '04-DESeq2_NormCounts'))
        .mix(pub(dge.dge_table, root, '05-DESeq2_DGE'))
    if (side_tables) {
        published = published
            .mix(pub(dge.contrasts, root, '05-DESeq2_DGE'))
            .mix(pub(dge.sample_table, root, '05-DESeq2_DGE'))
    }
    return published
}

def pub_vv(ch, root, name) {
    return ch.combine(root).combine(channel.value(name)).map { row ->
        def items = as_list(row)
        vv_log(items[-2], items[-1], items[0])
    }
}

def rewrite_figshare(url) {
    if ( url && url.toString().contains('figshare.com/ndownloader/files/') ) {
        def file_id = (url.toString() =~ /.*\/files\/([a-zA-Z0-9]+).*/)[0][1]
        return "https://api.figshare.com/v2/file/download/${file_id}"
    }
    return url
}

def is_remote_uri(p) {
    def s = p == null ? "" : p.toString().trim()
    return s.contains("://")
}

workflow RNASEQ {
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

    main:
        def ep = params.entry_point
        def microbes = params.mode == 'microbes'
        dge_script = "${projectDir}/bin/dge_deseq2.Rmd"
        ercc_notebook = "${projectDir}/bin/combined_ercc_analysis.ipynb"
        ch_multiqc_config = params.multiqc_config ? channel.fromPath( params.multiqc_config ) : channel.fromPath("NO_FILE")

        println "Entry point: '${ep}' mode: '${params.mode}'"
        if ( microbes && ep == 'genes_results' ) {
            error "entry_point 'genes_results' is eukaryotic-only (RSEM). Use counts_table for microbes."
        }

        STAGE(
            ch_outdir,
            dp_tools_plugin,
            accession,
            isa_archive_path,
            runsheet_path,
            api_url
        )
        ch_outdir = STAGE.out.ch_outdir
        samples = STAGE.out.samples
        samples_txt = STAGE.out.samples_txt
        runsheet_path = STAGE.out.runsheet_path
        isa_archive = STAGE.out.isa_archive
        osd_accession = STAGE.out.osd_accession
        ch_root = params.accession ? STAGE.out.glds_accession : channel.value('results')
        ch_published = STAGE.out.published

        if ( ep in ['counts_table', 'dge_table'] ) {
            samples | first | set { ch_meta }
        } else {
            samples | first | map { meta, _files -> meta } | set { ch_meta }
        }
        ch_meta | map { meta -> meta.organism_sci } | set { organism_sci }

        if ( is_remote_uri(params.reference_table) ) {
            FETCH_ANNOTATIONS_TABLE( annotations_csv_url_string )
            PARSE_ANNOTATIONS_TABLE( FETCH_ANNOTATIONS_TABLE.out.table, organism_sci )
        } else {
            PARSE_ANNOTATIONS_TABLE( annotations_csv_url_string.map { p -> file(p) }, organism_sci )
        }
        def ch_annot_src = params.gene_annotations_file
            ? channel.value(params.gene_annotations_file.toString())
            : PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url.map { url -> url ? url.toString() : "" }
        ch_annot_src.branch { p ->
            skip: !p
            remote: is_remote_uri(p)
            local: true
        }.set { ch_annot }
        FETCH_ANNOTATIONS( ch_annot.remote.map { rewrite_figshare(it) } )
        gene_annotations = FETCH_ANNOTATIONS.out.table
            .mix( ch_annot.local.map { p -> file(p) } )
            .mix( ch_annot.skip.map { file("${projectDir}/assets/empty") } )
            .first()

        if ( params.reference_fasta && params.reference_gtf ) {
            channel.value( params.reference_source ) | set { reference_source }
            channel.value( params.reference_version ) | set { reference_version }
            channel.value( params.reference_fasta ) | set { reference_fasta_url }
            channel.value( params.reference_gtf ) | set { reference_gtf_url }
        } else {
            reference_source = PARSE_ANNOTATIONS_TABLE.out.reference_source
            reference_version = PARSE_ANNOTATIONS_TABLE.out.reference_version
            reference_fasta_url = PARSE_ANNOTATIONS_TABLE.out.reference_fasta_url
            reference_gtf_url = PARSE_ANNOTATIONS_TABLE.out.reference_gtf_url
        }

        channel.empty() | set { genome_references }
        channel.empty() | set { genome_references_pre_ercc }
        channel.empty() | set { genome_bed }
        channel.empty() | set { ch_versions }

        if ( ep != 'dge_table' ) {
            reference_store_path
                .combine(organism_sci)
                .combine(reference_source)
                .combine(reference_version)
                .combine(reference_fasta_url)
                .combine(reference_gtf_url)
                .branch { store, org, src, ver, fasta, gtf ->
                    remote: is_remote_uri(fasta) || is_remote_uri(gtf)
                    local: true
                }
                .set { ch_ref }

            ch_ref.local.multiMap { store, org, src, ver, fasta, gtf ->
                store: store
                org: org
                src: src
                ver: ver
                fasta: file(fasta)
                gtf: file(gtf)
            }.set { ch_ref_local }
            COPY_REFERENCES(
                ch_ref_local.store,
                ch_ref_local.org,
                ch_ref_local.src,
                ch_ref_local.ver,
                ch_ref_local.fasta,
                ch_ref_local.gtf
            )
            ch_ref.remote.multiMap { store, org, src, ver, fasta, gtf ->
                store: store
                org: org
                src: src
                ver: ver
                fasta_url: fasta.toString()
                gtf_url: gtf.toString()
                fasta_local: is_remote_uri(fasta) ? file("${projectDir}/assets/empty") : file(fasta)
                gtf_local: is_remote_uri(gtf) ? file("${projectDir}/assets/empty") : file(gtf)
            }.set { ch_ref_remote }
            DOWNLOAD_REFERENCES(
                ch_ref_remote.store,
                ch_ref_remote.org,
                ch_ref_remote.src,
                ch_ref_remote.ver,
                ch_ref_remote.fasta_url,
                ch_ref_remote.gtf_url,
                ch_ref_remote.fasta_local,
                ch_ref_remote.gtf_local
            )
            genome_references_pre_subsample = COPY_REFERENCES.out.reference_files
                .mix(DOWNLOAD_REFERENCES.out.reference_files)

            if ( params.genome_subsample ) {
                SUBSAMPLE_GENOME( derived_store_path, organism_sci, genome_references_pre_subsample, reference_source, reference_version )
                SUBSAMPLE_GENOME.out.build | flatten | toList | set { genome_references_pre_ercc }
            } else {
                genome_references_pre_subsample | flatten | toList | set { genome_references_pre_ercc }
            }

            DOWNLOAD_ERCC( ch_meta.map { meta -> meta.has_ercc }, reference_store_path ).ifEmpty([file("ERCC92.fa"), file("ERCC92.gtf")]) | set { ch_maybe_ercc_refs }
            CONCAT_ERCC( reference_store_path, organism_sci, reference_source, reference_version, genome_references_pre_ercc, ch_maybe_ercc_refs, ch_meta.map { meta -> meta.has_ercc } )
                .ifEmpty { genome_references_pre_ercc.value } | first | set { genome_references }

            if ( ep in ['raw_reads', 'trimmed_reads', 'bam_files'] ) {
                GTF_TO_PRED(
                    derived_store_path,
                    organism_sci,
                    reference_source,
                    reference_version,
                    genome_references | map { refs -> refs[1] }
                )
                PRED_TO_BED(
                    derived_store_path,
                    organism_sci,
                    reference_source,
                    reference_version,
                    GTF_TO_PRED.out.genome_pred
                )
                genome_bed = PRED_TO_BED.out.genome_bed.first()
                ch_versions = ch_versions.mix(GTF_TO_PRED.out.versions).mix(PRED_TO_BED.out.versions)
            }
        }

        channel.empty() | set { raw_mqc_data }
        channel.empty() | set { raw_mqc_zip }
        channel.empty() | set { trimmed_mqc_data }
        channel.empty() | set { trimmed_mqc_zip }
        channel.empty() | set { sorted_bam }
        channel.empty() | set { bam_to_transcriptome }
        channel.empty() | set { reads_per_gene }
        channel.empty() | set { bam_only_files }
        channel.empty() | set { align_mqc_data }
        channel.empty() | set { align_mqc_zip }
        channel.empty() | set { genes_results }
        channel.empty() | set { rsem_publishables }
        channel.empty() | set { ercc_counts }
        channel.empty() | set { star_publishables }
        channel.empty() | set { count_mqc_data }
        channel.empty() | set { count_mqc_zip }
        channel.empty() | set { counts }
        channel.empty() | set { infer_mqc_data }
        channel.empty() | set { infer_mqc_zip }
        channel.empty() | set { genebody_mqc_data }
        channel.empty() | set { genebody_mqc_zip }
        channel.empty() | set { inner_mqc_data }
        channel.empty() | set { inner_mqc_zip }
        channel.empty() | set { readdist_mqc_data }
        channel.empty() | set { readdist_mqc_zip }
        channel.empty() | set { raw_fq_reports }
        channel.empty() | set { trimmed_fq_reports }
        channel.empty() | set { dge_table }
        channel.empty() | set { dge_table_rrnarm }
        channel.empty() | set { counts_rrnarm }
        channel.empty() | set { norm_counts }
        channel.empty() | set { qc_counts }

        nf_version = '"NEXTFLOW":\n    nextflow: '.concat("${nextflow.version}\n")
        ch_versions = channel.value(nf_version).mix(ch_versions)

        if ( ep in ['raw_reads', 'trimmed_reads'] ) {
            if ( ep == 'raw_reads' ) {
                RAW_FASTQC( STAGE.out.raw_reads )
                if ( !params.skip_vv ) {
                    VV_RAW_FASTQ( flatten_fastqs(STAGE.out.raw_reads) )
                    raw_fq_reports = VV_RAW_FASTQ.out.report
                }
                RAW_FASTQC.out.fastqc | map { qc -> [ qc[1], qc[2] ] }
                    | flatten
                    | collect
                    | set { raw_fastqc_zip }

                GET_MAX_READ_LENGTH( raw_fastqc_zip )
                max_read_length = GET_MAX_READ_LENGTH.out.length | map { n -> n.toString().toInteger() } | first

                TRIMGALORE( STAGE.out.raw_reads )
                trimmed_reads = TRIMGALORE.out.reads
                if ( !params.skip_vv ) {
                    VV_TRIMMED_FASTQ( flatten_fastqs(trimmed_reads) )
                    trimmed_fq_reports = VV_TRIMMED_FASTQ.out.report
                }

                TRIMMED_FASTQC( trimmed_reads )
                TRIMMED_FASTQC.out.fastqc | map { qc -> [ qc[1], qc[2] ] }
                    | flatten
                    | collect
                    | set { trimmed_fastqc_zip }

                RAW_READS_MULTIQC( samples_txt, raw_fastqc_zip, ch_multiqc_config, "raw_")
                TRIMMED_READS_MULTIQC( samples_txt, trimmed_fastqc_zip | concat( TRIMGALORE.out.reports ) | collect, ch_multiqc_config, "trimmed_")
                ch_published = ch_published
                    .mix( pub(RAW_FASTQC.out.fastqc, ch_root, '00-RawData/FastQC_Reports') )
                    .mix( pub(TRIMGALORE.out.reads, ch_root, '01-TG_Preproc/Fastq') )
                    .mix( pub(TRIMGALORE.out.reports, ch_root, '01-TG_Preproc/Trimming_Reports') )
                    .mix( pub(TRIMMED_FASTQC.out.fastqc, ch_root, '01-TG_Preproc/FastQC_Reports') )
                    .mix( pub_mqc(RAW_READS_MULTIQC.out, ch_root, '00-RawData/MultiQC_Reports') )
                    .mix( pub_mqc(TRIMMED_READS_MULTIQC.out, ch_root, '01-TG_Preproc/MultiQC_Reports') )

                raw_mqc_data = RAW_READS_MULTIQC.out.data
                raw_mqc_zip = RAW_READS_MULTIQC.out.zipped_data
                ch_versions = ch_versions.mix(RAW_FASTQC.out.versions).mix(TRIMGALORE.out.versions).mix(RAW_READS_MULTIQC.out.versions)
            } else {
                trimmed_reads = STAGE.out.trimmed_reads
                if ( !params.skip_vv ) {
                    VV_TRIMMED_FASTQ( flatten_fastqs(trimmed_reads) )
                    trimmed_fq_reports = VV_TRIMMED_FASTQ.out.report
                }

                TRIMMED_FASTQC( trimmed_reads )
                TRIMMED_FASTQC.out.fastqc | map { qc -> [ qc[1], qc[2] ] }
                    | flatten
                    | collect
                    | set { trimmed_fastqc_zip }

                GET_MAX_READ_LENGTH( trimmed_fastqc_zip )
                max_read_length = GET_MAX_READ_LENGTH.out.length | map { n -> n.toString().toInteger() } | first

                TRIMMED_READS_MULTIQC( samples_txt, trimmed_fastqc_zip, ch_multiqc_config, "trimmed_")
                ch_versions = ch_versions.mix(TRIMMED_READS_MULTIQC.out.versions)
                ch_published = ch_published
                    .mix( pub(TRIMMED_FASTQC.out.fastqc, ch_root, '01-TG_Preproc/FastQC_Reports') )
                    .mix( pub_mqc(TRIMMED_READS_MULTIQC.out, ch_root, '01-TG_Preproc/MultiQC_Reports') )
            }
            trimmed_mqc_data = TRIMMED_READS_MULTIQC.out.data
            trimmed_mqc_zip = TRIMMED_READS_MULTIQC.out.zipped_data

            if ( microbes ) {
                BUILD_BOWTIE2_INDEX( derived_store_path, organism_sci, reference_source, reference_version, genome_references, ch_meta )
                ALIGN_BOWTIE2( trimmed_reads, BUILD_BOWTIE2_INDEX.out.index_dir.first() )
                SAM_TO_BAM( ALIGN_BOWTIE2.out.sam )
                SORT_AND_INDEX_BAM( SAM_TO_BAM.out.bam )
                ALIGN_MULTIQC( samples_txt, ALIGN_BOWTIE2.out.alignment_logs.map { _m, f -> f } | collect, ch_multiqc_config, "align_")

                sorted_bam = SORT_AND_INDEX_BAM.out.sorted_bam
                bam_only_files = SORT_AND_INDEX_BAM.out.sorted_bam.map { row -> row[1] } | toSortedList()
                align_mqc_data = ALIGN_MULTIQC.out.data
                align_mqc_zip = ALIGN_MULTIQC.out.zipped_data
                ch_versions = ch_versions.mix(ALIGN_BOWTIE2.out.versions).mix(SAM_TO_BAM.out.versions)
                ch_published = ch_published
                    .mix( pub_by_meta(ALIGN_BOWTIE2.out.alignment_logs, ch_root, '02-Bowtie2_Alignment') )
                    .mix( pub_by_meta(ALIGN_BOWTIE2.out.unmapped_reads, ch_root, '02-Bowtie2_Alignment') )
                    .mix( pub_by_meta(SORT_AND_INDEX_BAM.out.sorted_bam, ch_root, '02-Bowtie2_Alignment') )
                    .mix( pub_mqc(ALIGN_MULTIQC.out, ch_root, '02-Bowtie2_Alignment/MultiQC_Reports') )
            } else {
                BUILD_STAR_INDEX( derived_store_path, organism_sci, reference_source, reference_version, genome_references, ch_meta, max_read_length )
                ALIGN_STAR( trimmed_reads, BUILD_STAR_INDEX.out.index_dir.first() )
                SORT_AND_INDEX_BAM( ALIGN_STAR.out.bam )
                ALIGN_MULTIQC( samples_txt, ALIGN_STAR.out.alignment_logs | collect, ch_multiqc_config, "align_")

                sorted_bam = SORT_AND_INDEX_BAM.out.sorted_bam
                bam_to_transcriptome = ALIGN_STAR.out.bam_to_transcriptome
                reads_per_gene = ALIGN_STAR.out.reads_per_gene
                bam_only_files = SORT_AND_INDEX_BAM.out.bam_only_files
                align_mqc_data = ALIGN_MULTIQC.out.data
                align_mqc_zip = ALIGN_MULTIQC.out.zipped_data
                ch_versions = ch_versions.mix(ALIGN_STAR.out.versions).mix(SORT_AND_INDEX_BAM.out.versions)
                ch_published = ch_published
                    .mix( ALIGN_STAR.out.publishables.combine(ch_root).flatMap { row ->
                        def items = as_list(row)
                        files_only(items[0..-2]).findAll { f -> f.name != 'versions.yml' }.collect { f ->
                            dest_file(items[-1], "02-STAR_Alignment/${f.parent.name}", f)
                        }
                    } )
                    .mix( pub_by_meta(SORT_AND_INDEX_BAM.out.sorted_bam, ch_root, '02-STAR_Alignment') )
                    .mix( pub_mqc(ALIGN_MULTIQC.out, ch_root, '02-STAR_Alignment/MultiQC_Reports') )
            }

            GENEBODY_COVERAGE( sorted_bam, genome_bed )
            INFER_EXPERIMENT( sorted_bam, genome_bed )
            INNER_DISTANCE( sorted_bam, genome_bed, max_read_length )
            READ_DISTRIBUTION( sorted_bam, genome_bed )

            infer_expt_out = INFER_EXPERIMENT.out.log | map { row -> row[1] } | collect
            ASSESS_STRANDEDNESS( infer_expt_out )
            if ( strandedness_set() ) {
                strandedness = channel.value(convert_strandedness(params.strandedness))
            } else {
                strandedness = ASSESS_STRANDEDNESS.out | map { f -> f.text.split(":")[0] } | first
            }

            INFER_EXPERIMENT_MULTIQC( samples_txt, INFER_EXPERIMENT.out.log | map { row -> row[1] } | collect, ch_multiqc_config, "infer_exp_")
            GENEBODY_COVERAGE_MULTIQC( samples_txt, GENEBODY_COVERAGE.out.log | map { row -> row[1] } | collect, ch_multiqc_config, "geneBody_cov_")
            INNER_DISTANCE_MULTIQC( samples_txt, INNER_DISTANCE.out.log | map { row -> row[1] } | collect, ch_multiqc_config, "inner_dist_")
            READ_DISTRIBUTION_MULTIQC( samples_txt, READ_DISTRIBUTION.out.log | map { row -> row[1] } | collect, ch_multiqc_config, "read_dist_")

            infer_mqc_data = INFER_EXPERIMENT_MULTIQC.out.data
            infer_mqc_zip = INFER_EXPERIMENT_MULTIQC.out.zipped_data
            genebody_mqc_data = GENEBODY_COVERAGE_MULTIQC.out.data
            genebody_mqc_zip = GENEBODY_COVERAGE_MULTIQC.out.zipped_data
            inner_mqc_data = INNER_DISTANCE_MULTIQC.out.data
            inner_mqc_zip = INNER_DISTANCE_MULTIQC.out.zipped_data
            readdist_mqc_data = READ_DISTRIBUTION_MULTIQC.out.data
            readdist_mqc_zip = READ_DISTRIBUTION_MULTIQC.out.zipped_data
            ch_versions = ch_versions.mix(INFER_EXPERIMENT.out.versions)
                .mix(GENEBODY_COVERAGE.out.versions)
                .mix(INNER_DISTANCE.out.versions)
                .mix(READ_DISTRIBUTION.out.versions)
            ch_published = ch_published
                .mix( pub_by_meta(GENEBODY_COVERAGE.out.all_output, ch_root, 'RSeQC_Analyses/02_geneBody_coverage') )
                .mix( pub(INFER_EXPERIMENT.out.log, ch_root, 'RSeQC_Analyses/03_infer_experiment') )
                .mix( pub_by_meta(INNER_DISTANCE.out.all_output, ch_root, 'RSeQC_Analyses/04_inner_distance') )
                .mix( pub(READ_DISTRIBUTION.out.log, ch_root, 'RSeQC_Analyses/05_read_distribution') )
                .mix( pub_mqc(INFER_EXPERIMENT_MULTIQC.out, ch_root, 'RSeQC_Analyses/MultiQC_Reports') )
                .mix( pub_mqc(GENEBODY_COVERAGE_MULTIQC.out, ch_root, 'RSeQC_Analyses/MultiQC_Reports') )
                .mix( pub_mqc(INNER_DISTANCE_MULTIQC.out, ch_root, 'RSeQC_Analyses/MultiQC_Reports') )
                .mix( pub_mqc(READ_DISTRIBUTION_MULTIQC.out, ch_root, 'RSeQC_Analyses/MultiQC_Reports') )
        } else if (ep == 'bam_files' && !strandedness_set()) {
            BAM_STRANDEDNESS(STAGE.out.bam_files, STAGE.out.genome_bam_files, genome_bed)
            strandedness = BAM_STRANDEDNESS.out.strandedness
            paired_counting_bams = BAM_STRANDEDNESS.out.counting_bams
            ch_versions = ch_versions.mix(BAM_STRANDEDNESS.out.versions)
            ch_published = ch_published.mix(pub(BAM_STRANDEDNESS.out.logs, ch_root, 'RSeQC_Analyses/03_infer_experiment'))
        } else {
            strandedness = channel.value(convert_strandedness(params.strandedness))
        }

        if ( !microbes && ep in ['raw_reads', 'trimmed_reads', 'bam_files', 'genes_results'] ) {
            if ( ep == 'genes_results' ) {
                genes_results = STAGE.out.genes_results
                rsem_counts = genes_results | map { row -> row[1] } | collect
                QUANTIFY_RSEM_GENES( samples_txt, rsem_counts )
            } else {
                if ( ep in ['raw_reads', 'trimmed_reads'] ) {
                    QUANTIFY_STAR_GENES( samples_txt, reads_per_gene | toSortedList, strandedness )
                    star_publishables = QUANTIFY_STAR_GENES.out.publishables
                    ch_published = ch_published.mix( pub(QUANTIFY_STAR_GENES.out.publishables, ch_root, '02-STAR_Alignment') )
                }

                def tx_bam = ep == 'bam_files' ? (strandedness_set() ? STAGE.out.bam_files : paired_counting_bams) : bam_to_transcriptome
                BUILD_RSEM_INDEX( derived_store_path, organism_sci, reference_source, reference_version, genome_references, ch_meta )
                COUNT_ALIGNED( tx_bam, BUILD_RSEM_INDEX.out.index_dir.first(), strandedness )
                genes_results = COUNT_ALIGNED.out.genes_results
                rsem_counts = COUNT_ALIGNED.out.counts | map { row -> row[1] } | collect
                QUANTIFY_RSEM_GENES( samples_txt, rsem_counts )

                COUNT_MULTIQC( samples_txt, rsem_counts, ch_multiqc_config, "RSEM_count_")
                count_mqc_data = COUNT_MULTIQC.out.data
                count_mqc_zip = COUNT_MULTIQC.out.zipped_data
                ch_versions = ch_versions.mix(COUNT_ALIGNED.out.versions).mix(COUNT_MULTIQC.out.versions)
                ch_published = ch_published
                    .mix( pub_by_meta(COUNT_ALIGNED.out.counts, ch_root, '03-RSEM_Counts') )
                    .mix( pub_mqc(COUNT_MULTIQC.out, ch_root, '03-RSEM_Counts/MultiQC_Reports') )
            }
            rsem_publishables = QUANTIFY_RSEM_GENES.out.publishables
            ch_published = ch_published.mix( pub(QUANTIFY_RSEM_GENES.out.publishables, ch_root, '03-RSEM_Counts') )
            qc_counts = QUANTIFY_RSEM_GENES.out.publishables
            ercc_counts = QUANTIFY_RSEM_GENES.out.publishables.map { unnorm, _nz -> unnorm }
        }

        if ( microbes && ep in ['raw_reads', 'trimmed_reads', 'bam_files'] ) {
            def bam_list = ep == 'bam_files' ? (strandedness_set() ? STAGE.out.bam_files : paired_counting_bams).map { row -> row[1] }.collect() : bam_only_files
            GET_GTF_FEATURES( genome_references )
            gtf_features = GET_GTF_FEATURES.out.gtf_features.map { f -> f.text.trim() }
            FEATURECOUNTS( ch_meta, genome_references, gtf_features, strandedness, bam_list )
            QUANTIFY_FEATURECOUNTS_GENES( samples_txt, FEATURECOUNTS.out.counts )
            qc_counts = QUANTIFY_FEATURECOUNTS_GENES.out.unnorm_counts
            COUNT_MULTIQC( samples_txt, FEATURECOUNTS.out.summary | collect, ch_multiqc_config, "FeatureCounts_")

            counts = FEATURECOUNTS.out.counts
            count_mqc_data = COUNT_MULTIQC.out.data
            count_mqc_zip = COUNT_MULTIQC.out.zipped_data
            ch_versions = ch_versions.mix(FEATURECOUNTS.out.versions).mix(COUNT_MULTIQC.out.versions)
            ch_published = ch_published
                .mix( pub(FEATURECOUNTS.out.counts, ch_root, '03-FeatureCounts') )
                .mix( pub(FEATURECOUNTS.out.summary, ch_root, '03-FeatureCounts') )
                .mix( pub(QUANTIFY_FEATURECOUNTS_GENES.out.num_non_zero_genes, ch_root, '03-FeatureCounts') )
                .mix( pub(QUANTIFY_FEATURECOUNTS_GENES.out.unnorm_counts, ch_root, '03-FeatureCounts') )
                .mix( pub_mqc(COUNT_MULTIQC.out, ch_root, '03-FeatureCounts/MultiQC_Reports') )
        }

        if ( !microbes && ep == 'counts_table' ) {
            ercc_counts = STAGE.out.counts_table
        }

        if ( !microbes && ep != 'dge_table' ) {
            ERCC_ANALYSIS(
                ch_meta.map { meta -> meta.has_ercc },
                STAGE.out.glds_accession.ifEmpty( params.accession ?: '' ),
                runsheet_path,
                channel.value( params.assay_suffix ),
                ercc_counts,
                ch_meta.map { meta -> meta.organism_sci },
                ercc_notebook
            )
            ch_published = ch_published.mix( pub(
                ERCC_ANALYSIS.out.notebook.mix(ERCC_ANALYSIS.out.html).mix(ERCC_ANALYSIS.out.error).mix(ERCC_ANALYSIS.out.tables),
                ch_root, 'ERCC_Analysis'
            ) )
        }

        dropped_unalignable = channel.value([])
        if ( ep == 'dge_table' ) {
            ANNOTATE_DGE_TABLE( gene_annotations, ch_meta, STAGE.out.dge_table )
            ch_versions = ch_versions.mix(ANNOTATE_DGE_TABLE.out.versions)
            ch_published = ch_published.mix( pub(ANNOTATE_DGE_TABLE.out.dge_table, ch_root, '05-DESeq2_DGE') )
        } else {
            EXTRACT_RRNA( organism_sci, genome_references | map { refs -> refs[1] } )
            if ( params.drop_unalignable.toString().toBoolean() && !microbes && ep in ['raw_reads', 'trimmed_reads', 'bam_files'] ) {
                FILTER_UNALIGNABLE( runsheet_path, count_mqc_data )
                dge_runsheet = FILTER_UNALIGNABLE.out.runsheet
                dropped_unalignable = FILTER_UNALIGNABLE.out.dropped
            } else {
                dge_runsheet = runsheet_path
            }
            SPLIT_ORGANISM_PART_RUNSHEETS( dge_runsheet )
            ch_dge_jobs = dge_jobs_from_manifest(
                SPLIT_ORGANISM_PART_RUNSHEETS.out.manifest,
                SPLIT_ORGANISM_PART_RUNSHEETS.out.runsheets
            )
            dge_rmd = channel.value( file(dge_script) )

            if ( ep == 'counts_table' ) {
                def counts_dir = microbes ? "/03-FeatureCounts" : "/03-RSEM_Counts"
                DGE_DESEQ2( dge_input(ch_dge_jobs, ch_meta, gene_annotations, STAGE.out.counts_table, false), dge_rmd )
                REMOVE_RRNA_COUNTS_TABLE( STAGE.out.counts_table, EXTRACT_RRNA.out.rrna_ids )
                counts_rrnarm = REMOVE_RRNA_COUNTS_TABLE.out.counts_rrnarm
                DGE_DESEQ2_RRNA_RM( dge_input(ch_dge_jobs, ch_meta, gene_annotations, counts_rrnarm, true), dge_rmd )
                ch_published = ch_published.mix( pub(REMOVE_RRNA_COUNTS_TABLE.out.counts_rrnarm, ch_root, counts_dir.replaceFirst('^/', '')) )
            } else if ( microbes ) {
                DGE_DESEQ2( dge_input(ch_dge_jobs, ch_meta, gene_annotations, counts, false), dge_rmd )
                REMOVE_RRNA_FEATURECOUNTS( counts, EXTRACT_RRNA.out.rrna_ids )
                counts_rrnarm = REMOVE_RRNA_FEATURECOUNTS.out.counts_rrnarm
                DGE_DESEQ2_RRNA_RM( dge_input(ch_dge_jobs, ch_meta, gene_annotations, counts_rrnarm, true), dge_rmd )
                REMOVE_RRNA_FC_TABLE(QUANTIFY_FEATURECOUNTS_GENES.out.unnorm_counts, EXTRACT_RRNA.out.rrna_ids)
                ch_published = ch_published
                    .mix( pub(REMOVE_RRNA_FEATURECOUNTS.out.counts_rrnarm, ch_root, '03-FeatureCounts') )
                    .mix( pub(REMOVE_RRNA_FC_TABLE.out.counts_rrnarm, ch_root, '03-FeatureCounts') )
            } else {
                DGE_DESEQ2( dge_input(ch_dge_jobs, ch_meta, gene_annotations, genes_results.map { row -> row[1] } | collect, false), dge_rmd )
                REMOVE_RRNA( EXTRACT_RRNA.out.rrna_ids, genes_results )
                counts_rrnarm = REMOVE_RRNA.out.genes_results_rrnarm.map { _m, f -> f } | toSortedList
                DGE_DESEQ2_RRNA_RM( dge_input(ch_dge_jobs, ch_meta, gene_annotations, counts_rrnarm, true), dge_rmd )
                REMOVE_RRNA_RSEM_TABLE(QUANTIFY_RSEM_GENES.out.publishables.map { unnorm, _nz -> unnorm }, EXTRACT_RRNA.out.rrna_ids)
                ch_published = ch_published
                    .mix( pub_by_meta(REMOVE_RRNA.out.genes_results_rrnarm, ch_root, '03-RSEM_Counts') )
                    .mix( pub(REMOVE_RRNA_RSEM_TABLE.out.counts_rrnarm, ch_root, '03-RSEM_Counts') )
            }

            dge_table = DGE_DESEQ2.out.dge_table
            dge_table_rrnarm = DGE_DESEQ2_RRNA_RM.out.dge_table
            norm_counts = DGE_DESEQ2.out.norm_counts
            ch_versions = ch_versions.mix(DGE_DESEQ2.out.versions)
            ch_published = ch_published
                .mix( pub_dge(DGE_DESEQ2.out, ch_root) )
                .mix( pub_dge(DGE_DESEQ2_RRNA_RM.out, ch_root, false) )
        }

        all_multiqc = raw_mqc_data
            | concat( trimmed_mqc_data )
            | concat( align_mqc_data )
            | concat( genebody_mqc_data )
            | concat( infer_mqc_data )
            | concat( inner_mqc_data )
            | concat( readdist_mqc_data )
            | concat( count_mqc_data )
            | collect

        if ( ep in ['raw_reads', 'trimmed_reads', 'bam_files', 'genes_results'] ) {
            PARSE_QC_METRICS(
                osd_accession,
                ch_meta,
                isa_archive.ifEmpty { [] },
                all_multiqc,
                qc_counts,
                runsheet_path
            )
            ch_published = ch_published.mix( pub(PARSE_QC_METRICS.out.file, ch_root, 'GeneLab') )
            ch_published = ch_published.mix( pub(PARSE_QC_METRICS.out.validation, ch_root, 'GeneLab') )
        }

        channel.empty() | set { vv_logs }

        if ( ep == 'raw_reads' ) {
            VV_RAW_READS( dp_tools_plugin, ch_outdir, ch_meta, runsheet_path, raw_mqc_zip, vv_fastq_reports(raw_fq_reports) )
            vv_logs = vv_logs.mix(VV_RAW_READS.out.log)
            ch_versions = ch_versions.mix(VV_RAW_READS.out.versions)
            ch_published = ch_published.mix( pub_vv(VV_RAW_READS.out.log, ch_root, 'VV_RAW_READS') )
        }
        if ( ep in ['raw_reads', 'trimmed_reads'] ) {
            VV_TRIMMED_READS( dp_tools_plugin, ch_outdir, ch_meta, runsheet_path, trimmed_mqc_zip, vv_fastq_reports(trimmed_fq_reports) )
            vv_logs = vv_logs.mix(VV_TRIMMED_READS.out.log)
            ch_versions = ch_versions.mix(VV_TRIMMED_READS.out.versions)
            ch_published = ch_published.mix( pub_vv(VV_TRIMMED_READS.out.log, ch_root, 'VV_TRIMMED_READS') )

            VV_RSEQC(
                dp_tools_plugin,
                ch_outdir,
                ch_meta,
                runsheet_path,
                genebody_mqc_zip,
                infer_mqc_zip,
                channel.empty() | mix(inner_mqc_zip) | collect | ifEmpty({ file("PLACEHOLDER") }),
                readdist_mqc_zip
            )
            vv_logs = vv_logs.mix(VV_RSEQC.out.log)
            ch_published = ch_published.mix( pub_vv(VV_RSEQC.out.log, ch_root, 'VV_RSEQC') )
        }
        if ( ep in ['raw_reads', 'trimmed_reads'] && !microbes ) {
            VV_STAR_ALIGNMENT(
                dp_tools_plugin,
                ch_outdir,
                ch_meta,
                runsheet_path,
                align_mqc_zip,
                star_publishables,
                bam_only_files | collect
            )
            vv_logs = vv_logs.mix(VV_STAR_ALIGNMENT.out.log)
            ch_published = ch_published.mix( pub_vv(VV_STAR_ALIGNMENT.out.log, ch_root, 'VV_ALIGNMENT') )
        }
        if ( ep in ['raw_reads', 'trimmed_reads'] && microbes ) {
            VV_BOWTIE2_ALIGNMENT(
                dp_tools_plugin,
                ch_outdir,
                ch_meta,
                runsheet_path,
                align_mqc_zip,
                bam_only_files | collect
            )
            vv_logs = vv_logs.mix(VV_BOWTIE2_ALIGNMENT.out.log)
            ch_published = ch_published.mix( pub_vv(VV_BOWTIE2_ALIGNMENT.out.log, ch_root, 'VV_ALIGNMENT') )
        }
        if ( ep in ['raw_reads', 'trimmed_reads', 'bam_files'] && !microbes ) {
            VV_RSEM_COUNTS(
                dp_tools_plugin,
                ch_outdir,
                ch_meta,
                runsheet_path,
                count_mqc_zip,
                counts_rrnarm,
                rsem_publishables.combine(REMOVE_RRNA_RSEM_TABLE.out.counts_rrnarm)
            )
            vv_logs = vv_logs.mix(VV_RSEM_COUNTS.out.log)
            ch_versions = ch_versions.mix(VV_RSEM_COUNTS.out.versions)
            ch_published = ch_published.mix( pub_vv(VV_RSEM_COUNTS.out.log, ch_root, 'VV_COUNTS') )
        }
        if ( ep in ['raw_reads', 'trimmed_reads', 'bam_files'] && microbes ) {
            VV_FEATURECOUNTS(
                dp_tools_plugin,
                ch_outdir,
                ch_meta,
                runsheet_path,
                count_mqc_zip,
                counts_rrnarm.combine(REMOVE_RRNA_FC_TABLE.out.counts_rrnarm)
            )
            vv_logs = vv_logs.mix(VV_FEATURECOUNTS.out.log)
            ch_versions = ch_versions.mix(VV_FEATURECOUNTS.out.versions)
            ch_published = ch_published.mix( pub_vv(VV_FEATURECOUNTS.out.log, ch_root, 'VV_COUNTS') )
        }
        if ( ep != 'dge_table' ) {
            VV_DGE_DESEQ2(
                dp_tools_plugin,
                ch_outdir,
                ch_meta,
                dge_runsheet,
                dge_table | collect,
                dge_table_rrnarm | collect,
                SPLIT_ORGANISM_PART_RUNSHEETS.out.stratify_by.map { f -> f.text.trim() }
            )
            vv_logs = vv_logs.mix(VV_DGE_DESEQ2.out.log)
            ch_versions = ch_versions.mix(VV_DGE_DESEQ2.out.versions)

            VV_CONCAT_FILTER( vv_logs | collect )
            ch_published = ch_published
                .mix( pub_vv(VV_DGE_DESEQ2.out.log, ch_root, 'VV_DESEQ2_ANALYSIS') )
                .mix( pub(VV_DGE_DESEQ2.out.log2fc_flags, ch_root, 'VV_Logs') )
                .mix( pub(VV_CONCAT_FILTER.out, ch_root, 'VV_Logs') )
        }

        ch_versions
            | unique
            | collectFile( newLine: true )
            | set { ch_final_software_versions }
        SOFTWARE_VERSIONS( ch_final_software_versions )
        ch_published = ch_published.mix( pub(SOFTWARE_VERSIONS.out.software_versions, ch_root, 'GeneLab') )

        if ( ep != 'dge_table' ) {
            GENERATE_PROTOCOL(
                ch_meta,
                strandedness,
                SOFTWARE_VERSIONS.out.software_versions_yaml,
                reference_source,
                reference_version,
                genome_references_pre_ercc,
                runsheet_path,
                dropped_unalignable
            )
            ch_published = ch_published.mix( pub(GENERATE_PROTOCOL.out.protocol, ch_root, 'GeneLab') )
        }

    emit:
        published = ch_published
}
