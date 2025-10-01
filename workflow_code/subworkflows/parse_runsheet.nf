// Adapted from Function: https://github.com/nf-core/rnaseq/blob/master/modules/local/process/samplesheet_check.nf
// Function to get list of [ meta, [ fastq_1_path, fastq_2_path ] ]
def get_runsheet_paths(LinkedHashMap row) {
    def meta = [:]
    def GENE_ID_TYPES = [
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

    meta.id = row["Sample Name"]
    meta.organism_sci = row.organism.replaceAll(" ","_").toLowerCase()
    meta.gene_id_type = GENE_ID_TYPES.get(meta.organism_sci, "gene_id")
    meta.paired_end = row.paired_end.toBoolean()
    meta.has_ercc = row.has_ERCC.toBoolean()

    // Extract factors
    meta.factors = row.findAll { key, value -> 
        key.startsWith("Factor Value[") && key.endsWith("]")
    }.collectEntries { key, value ->
        [(key[13..-2]): value] // Remove "Factor Value[" and "]"
    }

    def array = []
    def raw_reads = []
    raw_reads.add(file(row.read1_path))
    if (meta.paired_end) {
        raw_reads.add(file(row.read2_path))
      }
    array = [meta, raw_reads]
    return array
}

def mutate_to_single_end(it) {
    def new_meta = it[0].clone()  // Create a copy of the meta map
    new_meta.paired_end = false   // Set paired_end to false
    return [new_meta, [it[1][0]]] // Return only first read
}

// Truncate runsheet - Only used for debugging with params.limit_samples_to
// Only keep the first n sample rows
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
    
    main:
        sample_limit = params.limit_samples_to ? params.limit_samples_to : -1 // -1 in take means no limit

        // Only run truncation if there's an actual limit being applied
        if (sample_limit > 0) {
            // Create truncated runsheet using the process
            TRUNCATE_RUNSHEET(runsheet_path, sample_limit)
            ch_runsheet = TRUNCATE_RUNSHEET.out.truncated_runsheet
        } else {
            // Use the original runsheet directly
            ch_runsheet = runsheet_path
        }

        // Process samples from the runsheet
        ch_samples = ch_runsheet
            | splitCsv(header: true)
            | map { row -> get_runsheet_paths(row) }
            | map{ it -> params.force_single_end ? mutate_to_single_end(it) : it }

        ch_samples | set { ch_samples }

        // ch_samples | view

        // Validate consistency across samples
        ch_samples
            .map { meta, reads -> [meta.has_ercc, meta.paired_end, meta.organism_sci] }
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

        // Print autodetected processing metadata for the first sample
        ch_samples.take(1) | view { meta, reads -> 
            """Autodetected Processing Metadata:
            Has ERCC: ${meta.has_ercc}
            Paired End: ${meta.paired_end}
            Organism: ${meta.organism_sci}
            Gene ID Type: ${meta.gene_id_type}"""
        }
        // Check that all read file paths are unique
        ch_samples
            .flatMap { meta, reads -> reads }
            .collect()
            .map { all_reads ->
                def total_count = all_reads.size()
                def unique_count = all_reads.toSet().size()
                
                if (unique_count != total_count) {
                    throw new RuntimeException("ERROR: Duplicate read file paths detected. Please check the runsheet.")
                } else {
                    println "All ${unique_count} read file paths are unique."
                }
            }

    emit:
        samples = ch_samples
        runsheet = ch_runsheet
}


/* 
Entry point runsheet parsing workflows
*/

def get_trimmed_runsheet_paths(LinkedHashMap row) {
    def meta = [:]
    def GENE_ID_TYPES = [
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

    meta.id = row["Sample Name"]
    meta.organism_sci = row.organism.replaceAll(" ","_").toLowerCase()
    meta.gene_id_type = GENE_ID_TYPES.get(meta.organism_sci, "gene_id")
    meta.paired_end = row.paired_end.toBoolean()
    meta.has_ercc = row.has_ERCC.toBoolean()

    // Extract factors
    meta.factors = row.findAll { key, value -> 
        key.startsWith("Factor Value[") && key.endsWith("]")
    }.collectEntries { key, value ->
        [(key[13..-2]): value] // Remove "Factor Value[" and "]"
    }

    def array = []
    def trimmed_reads = []
    trimmed_reads.add(file(row.trimmed_read1_path))
    if (meta.paired_end) {
        trimmed_reads.add(file(row.trimmed_read2_path))
      }
    array = [meta, trimmed_reads]
    return array
}

workflow PARSE_TRIMMED_RUNSHEET {
    take:
        runsheet_path
    
    main:
        sample_limit = params.limit_samples_to ? params.limit_samples_to : -1 // -1 in take means no limit

        // Only run truncation if there's an actual limit being applied
        if (sample_limit > 0) {
            // Create truncated runsheet using the process
            TRUNCATE_RUNSHEET(runsheet_path, sample_limit)
            ch_runsheet = TRUNCATE_RUNSHEET.out.truncated_runsheet
        } else {
            // Use the original runsheet directly
            ch_runsheet = runsheet_path
        }

        // Process samples from the runsheet
        ch_samples = ch_runsheet
            | splitCsv(header: true)
            | map { row -> get_trimmed_runsheet_paths(row) }
            | map{ it -> params.force_single_end ? mutate_to_single_end(it) : it }

        ch_samples | set { ch_samples }

        // ch_samples | view

        // Validate consistency across samples
        ch_samples
            .map { meta, reads -> [meta.has_ercc, meta.paired_end, meta.organism_sci] }
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

        // Print autodetected processing metadata for the first sample
        ch_samples.take(1) | view { meta, reads -> 
            """Autodetected Processing Metadata:
            Has ERCC: ${meta.has_ercc}
            Paired End: ${meta.paired_end}
            Organism: ${meta.organism_sci}
            Gene ID Type: ${meta.gene_id_type}"""
        }

    emit:
        samples = ch_samples
        runsheet = ch_runsheet
}

def get_bam_runsheet_paths(LinkedHashMap row) {
    def meta = [:]
    def GENE_ID_TYPES = [
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

    meta.id = row["Sample Name"]
    meta.organism_sci = row.organism.replaceAll(" ","_").toLowerCase()
    meta.gene_id_type = GENE_ID_TYPES.get(meta.organism_sci, "gene_id")
    meta.paired_end = row.paired_end.toBoolean()
    meta.has_ercc = row.has_ERCC.toBoolean()

    // Extract factors
    meta.factors = row.findAll { key, value -> 
        key.startsWith("Factor Value[") && key.endsWith("]")
    }.collectEntries { key, value ->
        [(key[13..-2]): value] // Remove "Factor Value[" and "]"
    }

    def array = []
    def bam_files = []
    bam_files.add(file(row.bam_path))
    array = [meta, bam_files]
    return array
}

workflow PARSE_BAM_RUNSHEET {
    take:
        runsheet_path

    main:
        sample_limit = params.limit_samples_to ? params.limit_samples_to : -1 // -1 in take means no limit
  
        // Only run truncation if there's an actual limit being applied
        if (sample_limit > 0) {
            // Create truncated runsheet using the process
            TRUNCATE_RUNSHEET(runsheet_path, sample_limit)
            ch_runsheet = TRUNCATE_RUNSHEET.out.truncated_runsheet
        } else {
            // Use the original runsheet directly
            ch_runsheet = runsheet_path
        }

        // Process samples from the runsheet
        ch_samples = ch_runsheet
            | splitCsv(header: true)
            | map { row -> get_bam_runsheet_paths(row) }

        ch_samples | set { ch_samples }

        // ch_samples | view

        // Validate consistency across samples
        ch_samples
            .map { meta, reads -> [meta.has_ercc, meta.paired_end, meta.organism_sci] }
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

        // Print autodetected processing metadata for the first sample
        ch_samples.take(1) | view { meta, reads -> 
            """Autodetected Processing Metadata:
            Has ERCC: ${meta.has_ercc}
            Paired End: ${meta.paired_end}
            Organism: ${meta.organism_sci}
            Gene ID Type: ${meta.gene_id_type}"""
        }

    emit:
        samples = ch_samples
        runsheet = ch_runsheet
}

def get_genes_results_runsheet_paths(LinkedHashMap row) {
    def meta = [:]
    def GENE_ID_TYPES = [
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

    meta.id = row["Sample Name"]
    meta.organism_sci = row.organism.replaceAll(" ","_").toLowerCase()
    meta.gene_id_type = GENE_ID_TYPES.get(meta.organism_sci, "gene_id")
    meta.paired_end = row.paired_end.toBoolean()
    meta.has_ercc = row.has_ERCC.toBoolean()

    // Extract factors
    meta.factors = row.findAll { key, value -> 
        key.startsWith("Factor Value[") && key.endsWith("]")
    }.collectEntries { key, value ->
        [(key[13..-2]): value] // Remove "Factor Value[" and "]"
    }

    def array = []
    def genes_results_files = []
    genes_results_files.add(file(row.genes_results_path))
    array = [meta, genes_results_files]
    return array
}

workflow PARSE_GENES_RESULTS_RUNSHEET {
    take:
        runsheet_path

    main:
        sample_limit = params.limit_samples_to ? params.limit_samples_to : -1 // -1 in take means no limit
  
        // Only run truncation if there's an actual limit being applied
        if (sample_limit > 0) {
            // Create truncated runsheet using the process
            TRUNCATE_RUNSHEET(runsheet_path, sample_limit)
            ch_runsheet = TRUNCATE_RUNSHEET.out.truncated_runsheet
        } else {
            // Use the original runsheet directly
            ch_runsheet = runsheet_path
        }

        // Process samples from the runsheet
        ch_samples = ch_runsheet
            | splitCsv(header: true)
            | map { row -> get_genes_results_runsheet_paths(row) }

        ch_samples | set { ch_samples }

        // ch_samples | view

        // Validate consistency across samples
        ch_samples
            .map { meta, reads -> [meta.has_ercc, meta.paired_end, meta.organism_sci] }
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

        // Print autodetected processing metadata for the first sample
        ch_samples.take(1) | view { meta, reads -> 
            """Autodetected Processing Metadata:
            Has ERCC: ${meta.has_ercc}
            Paired End: ${meta.paired_end}
            Organism: ${meta.organism_sci}
            Gene ID Type: ${meta.gene_id_type}"""
        }

    emit:
        samples = ch_samples
        runsheet = ch_runsheet
}