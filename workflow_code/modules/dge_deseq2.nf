/*
 * DESeq2 Differential Gene Expression Analysis
 */
// ERCC counts are removed before normalization

process DGE_DESEQ2 {

    input:
        tuple val(meta), path(gene_annotations), val(output_label), path(runsheet_path), path(gene_counts)
        path("dge_deseq2.Rmd")

    output:
        tuple path("Normalized_Counts${output_label}${params.assay_suffix}.csv"),
              path(params.mode == "microbes" ? "FeatureCounts_Unnormalized_Counts${output_label}${params.assay_suffix}.csv" :
                   "RSEM_Unnormalized_Counts${output_label}${params.assay_suffix}.csv"), emit: norm_counts
        path("contrasts${output_label}${params.assay_suffix}.csv"), emit: contrasts, optional: true
        path("SampleTable${output_label}${params.assay_suffix}.csv"), emit: sample_table, optional: true
        path("differential_expression${output_label}${params.assay_suffix}.csv"), emit: dge_table
        path("VST_Counts${output_label}${params.assay_suffix}.csv"), emit: vst_norm_counts
        path("versions2.txt"), emit: versions

    script:
        def output_filename_label = output_label ?: ""
        def output_filename_suffix = params.assay_suffix ?: ""
        def microbes = params.mode == 'microbes' ? 'TRUE' : 'FALSE'
        if (!(params.dge_parallel_config in ['serial', 'multicore', 'snow']))
            error '--dge_parallel_config must be serial, multicore or snow'
        def cores = params.dge_parallel_cores
        if (!(cores.toString() ==~ /[1-9][0-9]*/) || !cores.toString().isInteger())
            error '--dge_parallel_cores must be a positive integer'
        def debug_dummy_counts = params.use_dummy_gene_counts ? 'TRUE'  : 'FALSE'

        // For counts_table entry point, pass the CSV file path; otherwise use directory/file path
        def input_counts_path = params.mode == 'microbes' ? gene_counts : "gene_counts"
        def use_counts_table = params.entry_point == "counts_table" ? 'TRUE' : 'FALSE'
        def counts_table_file = params.entry_point == "counts_table" ? gene_counts : ""

        """
        if [[ "${params.entry_point}" != "counts_table" ]] && [[ "${params.mode}" != "microbes" ]]; then
            mkdir -p gene_counts
            mv ${gene_counts} gene_counts/
        fi
        if [[ -s "${gene_annotations.toString()}" ]]; then
            ANNOT_ARG="${gene_annotations.toString()}"
        else
            ANNOT_ARG=""
        fi
        Rscript -e "rmarkdown::render('dge_deseq2.Rmd',
            output_file = 'DGE_DESeq2.html',
            output_dir = '\${PWD}',
            params = list(
                cpus = ${task.cpus},
                ${params.dge_parallel_config == 'serial' ? "" : "parallel_config = '${params.dge_parallel_config}',"}
                work_dir = '\${PWD}',
                output_directory = '\${PWD}',
                output_filename_label = '${output_filename_label}',
                output_filename_suffix = '${output_filename_suffix}',
                annotation_file_path = '\${ANNOT_ARG}',
                runsheet_path = '${runsheet_path}',
                microbes = ${microbes},
                gene_id_type = '${meta.gene_id_type}',
                input_counts = '${input_counts_path}',
                use_counts_table = ${use_counts_table},
                counts_table_file = '${counts_table_file}',
                DEBUG_MODE_LIMIT_GENES = FALSE,
                DEBUG_MODE_ADD_DUMMY_COUNTS = ${debug_dummy_counts},
                dge_filter_method = '${params.dge_filter_method}',
                dge_filter_sum_threshold = ${params.dge_filter_sum_threshold},
                dge_filter_sample_percent_threshold = ${params.dge_filter_sample_percent_threshold},
                dge_filter_sample_percent_max_samples = ${params.dge_filter_sample_percent_max_samples},
                dge_filter_min_samples_threshold = ${params.dge_filter_min_samples_threshold},
                dge_filter_count_per_sample_threshold = ${params.dge_filter_count_per_sample_threshold}
            ))"

        Rscript -e "versions <- c();
                    versions['R'] <- gsub(' .*', '', gsub('R version ', '', R.version\\\$version.string));
                    versions['BioConductor'] <- as.character(BiocManager::version());
                    pkg_list <- c('BiocParallel', 'DESeq2', 'tidyverse', 'dplyr', 'knitr', 'stringr', 'yaml');
                    if (${microbes} != TRUE) {
                        pkg_list <- c(pkg_list, 'tximport');
                    }
                    for(pkg in pkg_list) {
                        versions[pkg] <- as.character(packageVersion(pkg))
                    };
                    cat('"RNASEQ_DGE_DESEQ2":\\n',
                        paste0('    ', names(versions), ': ', versions, collapse='\\n'),
                        '\\n', sep='', file='versions2.txt')"
        """
}
