/*
    Downloads OSDR file list TSV for a dataset
*/

process GET_OSDR_FILE_LIST {
    
    input:
    val(osd_accession)

    output:
    path("*_file_list.tsv"), emit: file_list

    script:
    """
    python3 ${projectDir}/bin/osdr_downloader.py \\
        --osd ${osd_accession} \\
        --measurement "transcription profiling" \\
        --tech "RNA-Seq" \\
        --list \\
        --out .
    
    # Verify TSV was created
    if [ ! -f *_file_list.tsv ]; then
        echo "ERROR: Failed to download OSDR file list TSV"
        exit 1
    fi
    """
}

