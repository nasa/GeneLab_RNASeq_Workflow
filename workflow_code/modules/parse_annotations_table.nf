process PARSE_ANNOTATIONS_TABLE {

    input:
        path annotations_csv
        val organism_sci

    output:
        env 'fasta_url', emit: reference_fasta_url
        env 'gtf_url', emit: reference_gtf_url
        env 'gene_annotations_url', emit: gene_annotations_url
        env 'reference_source', emit: reference_source
        env 'reference_version', emit: reference_version

    script:
        def have_refs = (params.reference_fasta && params.reference_gtf) ? '1' : '0'
        def ep = params.entry_point ?: ''
        """
        parse_annotations_table.py \\
            '${annotations_csv}' \\
            '${organism_sci}' \\
            '${ep}' \\
            '${have_refs}'
        set -a
        . ./parsed.env
        set +a
        """
}
