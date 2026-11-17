process SPLIT_ORGANISM_PART_RUNSHEETS {
    input:
        path(runsheet)

    output:
        path("dge_jobs.tsv"), emit: manifest
        path("dge_part_*.csv"), emit: runsheets
        path("stratify_by.txt"), emit: stratify_by

    script:
        """
        split_organism_part_runsheets.py --runsheet ${runsheet}
        """
}
