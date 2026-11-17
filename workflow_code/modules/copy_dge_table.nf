process COPY_DGE_TABLE {

    input:
        path("?.csv")

    output:
        path("differential_expression${params.assay_suffix}.csv"), emit: dge_table

    script:
        def dge_table_name = "differential_expression${params.assay_suffix}.csv"
        """
        cp -P 1.csv ${dge_table_name}
        """
}

