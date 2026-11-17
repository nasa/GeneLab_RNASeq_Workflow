process COPY_REFERENCES {
  tag "Organism: ${organism_sci}, Reference Source: ${reference_source}${reference_source.toLowerCase().contains('ensembl') ? ', Reference Version: ' + reference_version : ''}"
  storeDir "${reference_store_path}/${reference_source}/${reference_source.toLowerCase().contains('ensembl') ? reference_version + '/' : ''}${organism_sci}"

  input:
    val(reference_store_path)
    val(organism_sci)
    val(reference_source)
    val(reference_version)
    path(fasta, stageAs: 'in_fasta/*')
    path(gtf, stageAs: 'in_gtf/*')

  output:
    tuple path("{*.fa,*.fna}"), path("*.gtf"), emit: reference_files

  script:
  """
  unpack_refs.sh "${fasta.toString()}" "${gtf.toString()}"
  """
}

process DOWNLOAD_REFERENCES {
  tag "Organism: ${organism_sci}, Reference Source: ${reference_source}${reference_source.toLowerCase().contains('ensembl') ? ', Reference Version: ' + reference_version : ''}"
  storeDir "${reference_store_path}/${reference_source}/${reference_source.toLowerCase().contains('ensembl') ? reference_version + '/' : ''}${organism_sci}"

  input:
    val(reference_store_path)
    val(organism_sci)
    val(reference_source)
    val(reference_version)
    val(fasta_url)
    val(gtf_url)
    path(fasta_local, stageAs: 'in_fasta/*')
    path(gtf_local, stageAs: 'in_gtf/*')

  output:
    tuple path("{*.fa,*.fna}"), path("*.gtf"), emit: reference_files

  script:
  """
  mkdir -p fetched
  if [[ "${fasta_url}" == *://* ]]; then
    fetch_uri.sh "${fasta_url}" "fetched/\$(basename "${fasta_url}")"
    fasta_src="fetched/\$(basename "${fasta_url}")"
  else
    fasta_src="${fasta_local.toString()}"
  fi
  if [[ "${gtf_url}" == *://* ]]; then
    fetch_uri.sh "${gtf_url}" "fetched/\$(basename "${gtf_url}")"
    gtf_src="fetched/\$(basename "${gtf_url}")"
  else
    gtf_src="${gtf_local.toString()}"
  fi
  unpack_refs.sh "\${fasta_src}" "\${gtf_src}"
  """
}
