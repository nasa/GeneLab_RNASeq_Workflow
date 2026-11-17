include { PUBLISH_METADATA } from '../modules/publish_metadata.nf'
include { PARSE_RUNSHEET } from './parse_runsheet.nf'
include { FETCH_ISA } from '../modules/fetch_isa.nf'
include { ISA_TO_RUNSHEET } from '../modules/isa_to_runsheet.nf'
include { GET_ACCESSIONS } from '../modules/get_accessions.nf'
include { STAGE_READS } from './stage_reads.nf'
include { GET_OSDR_FILE_LIST } from '../modules/get_osdr_file_list.nf'
include { UPDATE_OSDR_RUNSHEET } from '../modules/update_osdr_runsheet.nf'
include { COPY_BAMS as COPY_GENOME_BAMS } from '../modules/copy_bams.nf'
include { FETCH_BAM as FETCH_GENOME_BAM } from '../modules/fetch_remote.nf'
include { COPY_BAMS } from '../modules/copy_bams.nf'
include { COPY_GENES_RESULTS } from '../modules/copy_genes_results.nf'
include { COPY_COUNTS_TABLE } from '../modules/copy_counts_table.nf'
include { COPY_DGE_TABLE } from '../modules/copy_dge_table.nf'
include { FETCH_BAM; FETCH_GENES_RESULTS; FETCH_COUNTS_TABLE; FETCH_DGE_TABLE } from '../modules/fetch_remote.nf'
include { validateParameters } from 'plugin/nf-schema'

def is_remote_uri(p) {
    def s = p == null ? "" : p.toString().trim()
    return s.contains("://")
}

def entry_parse_types() {
    return [
        raw_reads:     'raw',
        trimmed_reads: 'trimmed',
        bam_files:     'bam',
        genes_results: 'genes',
        counts_table:  'counts',
        dge_table:     'dge'
    ]
}

def as_list(x) {
    if ( x == null ) return []
    if ( x instanceof Collection && !(x instanceof CharSequence) ) return x as List
    return [x]
}

def dest_file(root, dir, f) {
    return ["${root}/${dir}/${f.name}".toString(), f]
}

def files_only(items) {
    return items.collectMany { x -> x instanceof Map ? [] : as_list(x) }
}

def pub(ch, root, dir) {
    return ch.combine(root).combine(channel.value(dir)).flatMap { row ->
        def items = as_list(row)
        def d = items[-1]
        def r = items[-2]
        files_only(items[0..-3]).collect { f -> dest_file(r, d, f) }
    }
}

def samples_txt_from(ch_pairs) {
    return ch_pairs
        .map { pair -> pair[0].id }
        .collectFile(name: "samples.txt", sort: true, newLine: true)
}

def remote_files(ch) {
    ch.branch { _meta, files ->
        remote: files && files[0]?.toString()?.contains("://")
        local: true
    }
}

def remote_path(ch) {
    ch.branch { p ->
        remote: is_remote_uri(p)
        local: true
    }
}

def pub_sample(ch, root, dir) {
    return ch.combine(root).flatMap { row ->
        def items = as_list(row)
        def meta = items[0]
        files_only(items[1..-2]).collect { f -> dest_file(items[-1], "${dir}/${meta.id}", f) }
    }
}

workflow STAGE {
    take:
        ch_outdir
        dp_tools_plugin
        accession
        isa_archive_path
        runsheet_path
        api_url

    main:
        def ep = params.entry_point
        def parse_type = entry_parse_types()[ep]

        channel.empty() | set { osd_accession }
        channel.empty() | set { glds_accession }

        if ( accession ) {
            GET_ACCESSIONS( accession, api_url )
            osd_accession = GET_ACCESSIONS.out.accessions_txt.map { txt -> txt.readLines()[0].trim() }
            glds_accession = GET_ACCESSIONS.out.accessions_txt.map { txt -> txt.readLines()[1].trim() }
            ch_outdir = ch_outdir.combine(glds_accession).map { outdir, glds -> "$outdir/$glds" }
        } else {
            ch_outdir = ch_outdir.map { dir -> dir + "/results" }
        }
        ch_outdir = ch_outdir.first()
        ch_root = accession ? glds_accession : channel.value('results')

        channel.empty() | set { published }
        channel.empty() | set { isa_archive }
        if ( runsheet_path == null ) {
            if ( isa_archive_path == null ) {
                FETCH_ISA( osd_accession, glds_accession )
                isa_archive = FETCH_ISA.out.isa_archive
            } else {
                isa_archive = isa_archive_path
            }
            ISA_TO_RUNSHEET( osd_accession, glds_accession, isa_archive, dp_tools_plugin )
            runsheet_path = ISA_TO_RUNSHEET.out.runsheet
        } else if ( isa_archive_path != null ) {
            isa_archive = isa_archive_path
        }

        if ( params.validate_params ) {
            validateParameters()
        }

        channel.empty() | set { samples }
        channel.empty() | set { samples_txt }
        channel.empty() | set { raw_reads }
        channel.empty() | set { trimmed_reads }
        channel.empty() | set { bam_files }
        channel.empty() | set { genome_bam_files }
        channel.empty() | set { genes_results }
        channel.empty() | set { counts_table }
        channel.empty() | set { dge_table }

        def no_runsheet = (params.runsheet == null)
        def counts_override = (ep == 'counts_table' && params.counts_table_path)
        def dge_override = (ep == 'dge_table' && params.dge_table_path)
        def update_osdr = no_runsheet && ep != 'raw_reads' && !counts_override && !dge_override

        if ( update_osdr ) {
            GET_OSDR_FILE_LIST( osd_accession )
            UPDATE_OSDR_RUNSHEET(
                runsheet_path,
                GET_OSDR_FILE_LIST.out.file_list,
                glds_accession,
                ep
            )
            runsheet_path = UPDATE_OSDR_RUNSHEET.out.runsheet
        }

        def type = (counts_override || dge_override) ? 'meta' : parse_type
        PARSE_RUNSHEET( runsheet_path, channel.value(type) )
        runsheet_path = PARSE_RUNSHEET.out.runsheet

        if ( counts_override ) {
            samples = PARSE_RUNSHEET.out.samples.map { meta, _files -> meta }
            if ( is_remote_uri(params.counts_table_path) ) {
                FETCH_COUNTS_TABLE( channel.value(params.counts_table_path.toString()) )
                counts_table = FETCH_COUNTS_TABLE.out.counts_table
            } else {
                COPY_COUNTS_TABLE( file(params.counts_table_path) )
                counts_table = COPY_COUNTS_TABLE.out.counts_table
            }
            published = published.mix( pub(counts_table, ch_root, params.mode == 'microbes' ? '03-FeatureCounts' : '03-RSEM_Counts') )
        } else if ( dge_override ) {
            samples = PARSE_RUNSHEET.out.samples.map { meta, _files -> meta }
            if ( is_remote_uri(params.dge_table_path) ) {
                FETCH_DGE_TABLE( channel.value(params.dge_table_path.toString()) )
                dge_table = FETCH_DGE_TABLE.out.dge_table
            } else {
                COPY_DGE_TABLE( file(params.dge_table_path) )
                dge_table = COPY_DGE_TABLE.out.dge_table
            }
            published = published.mix( pub(dge_table, ch_root, '05-DESeq2_DGE') )
        } else if ( ep == 'raw_reads' ) {
            samples = PARSE_RUNSHEET.out.samples
            STAGE_READS( samples, channel.value("raw") )
            raw_reads = STAGE_READS.out.reads
            samples_txt = STAGE_READS.out.samples_txt
            published = published.mix( pub(STAGE_READS.out.reads, ch_root, '00-RawData/Fastq') )
        } else if ( ep == 'trimmed_reads' ) {
            samples = PARSE_RUNSHEET.out.samples
            STAGE_READS( samples, channel.value("trimmed") )
            trimmed_reads = STAGE_READS.out.reads
            samples_txt = STAGE_READS.out.samples_txt
            published = published.mix( pub(STAGE_READS.out.reads, ch_root, '01-TG_Preproc/Fastq') )
        } else if ( ep == 'bam_files' ) {
            samples = PARSE_RUNSHEET.out.samples
            ch_bam = remote_files(samples)
            FETCH_BAM( ch_bam.remote )
            COPY_BAMS( ch_bam.local )
            bam_files = FETCH_BAM.out.bam_files.mix( COPY_BAMS.out.bam_files )
            if (params.strandedness == 'auto' && params.mode != 'microbes') {
                ch_genome = remote_files(PARSE_RUNSHEET.out.genome_samples)
                FETCH_GENOME_BAM(ch_genome.remote)
                COPY_GENOME_BAMS(ch_genome.local)
                genome_bam_files = FETCH_GENOME_BAM.out.bam_files.mix(COPY_GENOME_BAMS.out.bam_files)
                published = published.mix(pub_sample(genome_bam_files, ch_root, '02-STAR_Alignment'))
            }
            samples_txt = samples_txt_from(bam_files)
            published = published.mix( pub_sample(bam_files, ch_root, params.mode == 'microbes' ? '02-Bowtie2_Alignment' : '02-STAR_Alignment') )
        } else if ( ep == 'genes_results' ) {
            samples = PARSE_RUNSHEET.out.samples
            ch_genes = remote_files(samples)
            FETCH_GENES_RESULTS( ch_genes.remote )
            COPY_GENES_RESULTS( ch_genes.local )
            genes_results = FETCH_GENES_RESULTS.out.genes_results.mix( COPY_GENES_RESULTS.out.genes_results )
            samples_txt = samples_txt_from(genes_results)
            published = published.mix( pub_sample(genes_results, ch_root, '03-RSEM_Counts') )
        } else if ( ep == 'counts_table' ) {
            samples = PARSE_RUNSHEET.out.samples.map { meta, _files -> meta }
            ch_counts = remote_path(PARSE_RUNSHEET.out.table)
            FETCH_COUNTS_TABLE(ch_counts.remote)
            COPY_COUNTS_TABLE(ch_counts.local.map { p -> file(p) })
            counts_table = FETCH_COUNTS_TABLE.out.counts_table.mix( COPY_COUNTS_TABLE.out.counts_table )
            published = published.mix( pub(counts_table, ch_root, params.mode == 'microbes' ? '03-FeatureCounts' : '03-RSEM_Counts') )
        } else if ( ep == 'dge_table' ) {
            samples = PARSE_RUNSHEET.out.samples.map { meta, _files -> meta }
            ch_dge = remote_path(PARSE_RUNSHEET.out.table)
            FETCH_DGE_TABLE(ch_dge.remote)
            COPY_DGE_TABLE(ch_dge.local.map { p -> file(p) })
            dge_table = FETCH_DGE_TABLE.out.dge_table.mix( COPY_DGE_TABLE.out.dge_table )
            published = published.mix( pub(dge_table, ch_root, '05-DESeq2_DGE') )
        }

        PUBLISH_METADATA( runsheet_path, isa_archive.ifEmpty { [] } )
        published = published.mix( pub(PUBLISH_METADATA.out, ch_root, 'Metadata') )

    emit:
        ch_outdir       = ch_outdir
        samples         = samples
        samples_txt     = samples_txt
        runsheet_path   = runsheet_path
        isa_archive     = isa_archive
        osd_accession   = osd_accession
        glds_accession  = glds_accession
        raw_reads       = raw_reads
        trimmed_reads   = trimmed_reads
        genome_bam_files = genome_bam_files
        bam_files       = bam_files
        genes_results   = genes_results
        counts_table    = counts_table
        dge_table       = dge_table
        published       = published
}
