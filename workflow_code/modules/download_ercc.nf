process DOWNLOAD_ERCC {
  storeDir "${reference_store_path}/ERCC_thermofisher"

  input:
    val(has_ercc)
    val(reference_store_path)

  output:
    tuple path("ERCC92.fa"), path("ERCC92.gtf")

  when:
    has_ercc

  script:
    """
    fetch_uri.sh --insecure \
        https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip \
        ERCC92.zip
    unzip ERCC92.zip
    """
}