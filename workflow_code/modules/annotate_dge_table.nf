
process ANNOTATE_DGE_TABLE {

    publishDir "${ publishdir }/05-DESeq2_DGE",
        pattern: "differential_expression${params.assay_suffix}.csv", 
        mode: params.publish_dir_mode

    input:
        val(publishdir)
        val(gene_annotations_url)
        val(meta)
        path("?.csv")
        path("annotate_dge_table.Rmd")

    output:
        path("differential_expression${params.assay_suffix}.csv"),        emit: dge_table
        path("versions.txt"),                                             emit: versions

    script:
        def output_filename_suffix = params.assay_suffix ?: ""

        """
        Rscript -e "rmarkdown::render('annotate_dge_table.Rmd', 
            output_file = 'Annotate_DGE_Table.html',
            output_dir = '\${PWD}',
            params = list(
                work_dir = '\${PWD}',
                output_directory = '\${PWD}',
                output_filename_suffix = '${output_filename_suffix}',
                annotation_file_path = '${gene_annotations_url}',
                gene_id_type = '${meta.gene_id_type}',
                input_table_path = '1.csv'
            ))"

        Rscript -e "versions <- c(); 
                    versions['R'] <- gsub(' .*', '', gsub('R version ', '', R.version\\\$version.string));
                    pkg_list <- c('dplyr', 'knitr');
                    for(pkg in pkg_list) {
                        versions[pkg] <- as.character(packageVersion(pkg))
                    };
                    cat('"ANNOTATE_DGE_TABLE":\\n', 
                        paste0('    ', names(versions), ': ', versions, collapse='\\n'), 
                        '\\n', sep='', file='versions.txt')"
        """
}