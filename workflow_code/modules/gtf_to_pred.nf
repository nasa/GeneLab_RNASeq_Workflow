process GTF_TO_PRED {
    // Converts reference gtf into pred 
    storeDir "${ derived_store_path }/Genome_GTF_BED_Files/${reference_source}/${reference_source.toLowerCase().contains('ensembl') ? reference_version + '/' : ''}${organism_sci}"
    
    input:
        val(derived_store_path)
        val(organism_sci)
        val(reference_source) // Used for defining storage location 
        val(reference_version) // Used for defining storage location 
        path(genome_gtf)

    output:
        path("${ genome_gtf }.genePred"), emit: genome_pred
        path("versions_gtfToGenePred.yml"), emit: versions

    script:
    """
    gtfToGenePred -geneNameAsName2 -ignoreGroupsWithoutExons ${ genome_gtf } ${ genome_gtf }.genePred

    echo '"${task.process}":' > versions_gtfToGenePred.yml
    echo "    gtfToGenePred: 469" >> versions_gtfToGenePred.yml
    """
}