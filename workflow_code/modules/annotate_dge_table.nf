
process ANNOTATE_DGE_TABLE {

    publishDir "${ publishdir }/05-DESeq2_DGE",
        pattern: "differential_expression${params.assay_suffix}.csv", 
        mode: params.publish_dir_mode

    input:
        val(publishdir)
        val(gene_annotations_url)
        val(meta)
        path("?.csv")

    output:
        path("differential_expression${params.assay_suffix}.csv"),        emit: dge_table
        path("versions.txt"),                                             emit: versions

    script:
        def output_filename_suffix = params.assay_suffix ?: ""

        """
        annotate_dge_table.R \
            '1.csv' \
            '${gene_annotations_url}' \
            '${meta.gene_id_type}' \
            'differential_expression${output_filename_suffix}.csv'

        Rscript -e "versions <- c(); 
                    versions['R'] <- gsub(' .*', '', gsub('R version ', '', R.version\\\$version.string));
                    pkg_list <- c('dplyr');
                    for(pkg in pkg_list) {
                        versions[pkg] <- as.character(packageVersion(pkg))
                    };
                    cat('"ANNOTATE_DGE_TABLE":\\n', 
                        paste0('    ', names(versions), ': ', versions, collapse='\\n'), 
                        '\\n', sep='', file='versions.txt')"
        """
}