// Adapted from Function: https://github.com/nf-core/rnaseq/blob/master/modules/local/process/samplesheet_check.nf
// Function to get list of [ meta, [ fastq_1_path, fastq_2_path ] ]
def gene_id_types() {
    return [
        // Mammals
        "homo_sapiens": "ENSEMBL",
        "mus_musculus": "ENSEMBL",
        "rattus_norvegicus": "ENSEMBL",

        // Other Vertebrates
        "danio_rerio": "ENSEMBL",
        "oryzias_latipes": "ENSEMBL",

        // Invertebrates
        "caenorhabditis_elegans": "ENSEMBL",
        "drosophila_melanogaster": "ENSEMBL",

        // Plants
        "arabidopsis_thaliana": "TAIR",
        "brachypodium_distachyon": "ENSEMBL",
        "oryza_sativa": "ENSEMBL",

        // Microbes
        "bacillus_subtilis": "ENSEMBL",
        "escherichia_coli": "ENSEMBL",
        "lactobacillus_acidophilus": "LOCUS",
        "mycobacterium_marinum": "LOCUS",
        "pseudomonas_aeruginosa": "LOCUS",
        "salmonella_enterica": "ENSEMBL",
        "saccharomyces_cerevisiae": "ENSEMBL",
        "serratia_liquefaciens": "LOCUS",
        "staphylococcus_aureus": "LOCUS",
        "streptococcus_mutans": "LOCUS",
        "vibrio_fischeri": "LOCUS"
    ]
}

// type -> runsheet columns. paired extras only used when meta.paired_end.
// counts/dge are dataset-wide paths (not per-sample files).
def runsheet_cols() {
    return [
        raw:      [files: ['read1_path'],              paired: ['read2_path']],
        trimmed:  [files: ['trimmed_read1_path'],      paired: ['trimmed_read2_path']],
        bam:      [files: ['bam_path']],
        genome_bam: [files: ['genome_bam_path']],
        genes:    [files: ['genes_results_path']],
        counts:   [table: 'counts_table_path'],
        dge:      [table: 'dge_table_path'],
        meta:     [:]
    ]
}

def row_to_meta(LinkedHashMap row) {
    def meta = [:]
    meta.id = row["Sample Name"]
    meta.organism_sci = row.organism.replaceAll(" ","_").toLowerCase()
    meta.gene_id_type = gene_id_types().get(meta.organism_sci, "gene_id")
    meta.paired_end = row.paired_end.toBoolean()
    meta.has_ercc = row.has_ERCC.toBoolean()
    meta.factors = row.findAll { key, _value ->
        key.startsWith("Factor Value[") && key.endsWith("]")
    }.collectEntries { key, value ->
        [(key[13..-2]): value]
    }
    return meta
}

def is_remote_uri(p) {
    def s = p == null ? "" : p.toString().trim()
    return s.contains("://")
}

def stage_runsheet_path(p) {
    def s = p.toString().trim()
    return is_remote_uri(s) ? s : file(s)
}

def require_col(LinkedHashMap row, String col, String type) {
    def v = row.containsKey(col) ? row[col] : null
    if (v == null || v.toString().trim() == "") {
        throw new RuntimeException(
            "Sample '${row['Sample Name']}': --entry_point ${params.entry_point} (type=${type}) needs column '${col}'"
        )
    }
    return v
}

def get_runsheet_paths(LinkedHashMap row, String type) {
    def meta = row_to_meta(row)
    if (type == 'genome_bam') meta.genome_bam = true
    def cols = runsheet_cols()
    if (!cols.containsKey(type)) {
        throw new RuntimeException("Unknown runsheet type '${type}'. Expected one of: ${cols.keySet()}")
    }
    def spec = cols[type]
    // meta: ISA / --counts_table_path override — sample metadata only
    if (!spec.files && !spec.table) {
        return [meta, []]
    }
    if (spec.table) {
        require_col(row, spec.table, type)
        return [meta, []]
    }
    def files = []
    spec.files.each { col -> files.add(stage_runsheet_path(require_col(row, col, type))) }
    if (meta.paired_end && spec.paired) {
        spec.paired.each { col -> files.add(stage_runsheet_path(require_col(row, col, type))) }
    }
    def remote = files.collect { f -> is_remote_uri(f) }
    if ( remote.any { r -> r } && remote.any { r -> !r } ) {
        throw new RuntimeException("Sample '${meta.id}': mix of remote and local paths in the runsheet")
    }
    return [meta, files]
}

def table_path_from_row(LinkedHashMap row, String type) {
    def spec = runsheet_cols()[type]
    spec?.table ? require_col(row, spec.table, type) : null
}

def mutate_to_single_end(sample) {
    def new_meta = sample[0].clone()
    new_meta.paired_end = false
    return [new_meta, [sample[1][0]]]
}

process TRUNCATE_RUNSHEET {

    input:
    path(runsheet)
    val(limit)

    output:
    path "${runsheet.baseName}_truncated.csv", emit: truncated_runsheet

    script:
    """
    head -n 1 ${runsheet} > ${runsheet.baseName}_truncated.csv
    if [ ${limit} -gt 0 ]; then
        tail -n +2 ${runsheet} | head -n ${limit} >> ${runsheet.baseName}_truncated.csv
    else
        tail -n +2 ${runsheet} >> ${runsheet.baseName}_truncated.csv
    fi
    """
}

workflow PARSE_RUNSHEET {
    take:
        runsheet_path
        type

    main:
        sample_limit = params.limit_samples_to ? params.limit_samples_to : -1

        if (sample_limit > 0) {
            TRUNCATE_RUNSHEET(runsheet_path, sample_limit)
            ch_runsheet = TRUNCATE_RUNSHEET.out.truncated_runsheet
        } else {
            ch_runsheet = runsheet_path
        }

        ch_rows = ch_runsheet
            | splitCsv(header: true, quote: '"')
            | combine(type)

        ch_samples = ch_rows
            | map { row, parse_type ->
                def item = get_runsheet_paths(row, parse_type)
                (params.force_single_end && (parse_type == 'raw' || parse_type == 'trimmed')) ? mutate_to_single_end(item) : item
            }

        ch_samples
            .map { meta, _files -> [meta.has_ercc, meta.paired_end, meta.organism_sci] }
            .unique()
            .count()
            .subscribe { count ->
                if (count > 1) {
                    log.error "ERROR: Inconsistent metadata across samples. Please check the runsheet."
                    exit 1
                } else {
                    println "Metadata consistency check passed."
                }
            }

        ch_samples.take(1) | view { meta, _files ->
            """Autodetected Processing Metadata:
            Has ERCC: ${meta.has_ercc}
            Paired End: ${meta.paired_end}
            Organism: ${meta.organism_sci}
            Gene ID Type: ${meta.gene_id_type}"""
        }

        ch_samples
            .flatMap { _meta, files -> files }
            .collect()
            .map { all_files ->
                if (all_files.size() == 0) {
                    return
                }
                def unique_count = all_files.toSet().size()
                if (unique_count != all_files.size()) {
                    throw new RuntimeException("ERROR: Duplicate file paths detected. Please check the runsheet.")
                } else {
                    println "All ${unique_count} file paths are unique."
                }
            }

        table = ch_rows
            | map { row, parse_type -> table_path_from_row(row, parse_type) }
            | filter { path -> path != null && path != "" }
            | unique

        genome_samples = ch_rows
            .filter { row, parse_type -> parse_type == 'bam' && params.mode != 'microbes' && params.strandedness == 'auto' }
            .map { row, _type -> get_runsheet_paths(row, 'genome_bam') }

        genome_samples.flatMap { meta, files -> files }.collect().map { paths ->
            if (paths.size() != paths.toSet().size()) error "BAM auto strandedness requires a distinct genome_bam_path for each sample"
        }

        if (params.entry_point == 'bam_files' && params.strandedness == 'auto') {
            ch_samples.map { meta, _files -> meta.id }.collect().map { ids ->
                if (ids.size() != ids.toSet().size()) error "BAM auto strandedness requires unique Sample Name values"
            }
        }

    emit:
        genome_samples = genome_samples
        samples = ch_samples
        table = table
        runsheet = ch_runsheet
}
