process FETCH_READS {
    tag "Sample: ${meta.id}"

    input:
        tuple val(meta), val(uris)
        val(type)

    output:
        tuple val(meta), path("${meta.id}${params.assay_suffix}*.gz"), emit: reads

    script:
        def suffix = type == "raw" ? "raw" : "trimmed"
        def out1 = "${meta.id}${params.assay_suffix}_R1_${suffix}.fastq.gz"
        def out2 = "${meta.id}${params.assay_suffix}_R2_${suffix}.fastq.gz"
        def out_se = "${meta.id}${params.assay_suffix}_${suffix}.fastq.gz"
        def u2 = (uris.size() > 1) ? uris[1] : ""
        def nrec = params.truncate_to ? params.truncate_to.toInteger() : 0
        if ( meta.paired_end && nrec ) {
        """
        fastq_first_n.sh ${nrec} ${out1} '${uris[0]}'
        fastq_first_n.sh ${nrec} ${out2} '${u2}'
        """
        } else if ( meta.paired_end ) {
        """
        fetch_uri.sh '${uris[0]}' ${out1}
        fetch_uri.sh '${u2}' ${out2}
        """
        } else if ( nrec ) {
        """
        fastq_first_n.sh ${nrec} ${out_se} '${uris[0]}'
        """
        } else {
        """
        fetch_uri.sh '${uris[0]}' ${out_se}
        """
        }
}

process FETCH_BAM {
    tag "Sample: ${meta.id}"

    input:
        tuple val(meta), val(uris)

    output:
        tuple val(meta), path("${meta.id}*.bam"), emit: bam_files

    script:
        def bam_name = meta.genome_bam ? "${meta.id}${params.assay_suffix}_Aligned.out.bam" : params.mode == "microbes" ?
            "${meta.id}${params.assay_suffix}.bam" :
            "${meta.id}${params.assay_suffix}_Aligned.toTranscriptome.out.bam"
        """
        fetch_uri.sh '${uris[0]}' ${bam_name}
        """
}

process FETCH_GENES_RESULTS {
    tag "Sample: ${meta.id}"

    input:
        tuple val(meta), val(uris)

    output:
        tuple val(meta), path("${meta.id}${params.assay_suffix}.genes.results"), emit: genes_results

    script:
        """
        fetch_uri.sh '${uris[0]}' ${meta.id}${params.assay_suffix}.genes.results
        """
}

process FETCH_TABLE {
    input:
        val(uri)

    output:
        path("1.csv"), emit: table

    script:
        """
        fetch_uri.sh '${uri}' 1.csv
        """
}

process FETCH_COUNTS_TABLE {
    input:
        val(uri)

    output:
        path("*Unnormalized_Counts*.csv"), emit: counts_table

    script:
        def counts_name = params.mode == "microbes" ?
            "FeatureCounts_Unnormalized_Counts${params.assay_suffix}.csv" :
            "RSEM_Unnormalized_Counts${params.assay_suffix}.csv"
        """
        fetch_uri.sh '${uri}' ${counts_name}
        """
}

process FETCH_DGE_TABLE {
    input:
        val(uri)

    output:
        path("differential_expression${params.assay_suffix}.csv"), emit: dge_table

    script:
        """
        fetch_uri.sh '${uri}' differential_expression${params.assay_suffix}.csv
        """
}
