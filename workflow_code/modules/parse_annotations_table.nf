process PARSE_ANNOTATIONS_TABLE {
  // Extracts data from this kind of table: 
  // https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv

  input:
    val(annotations_csv_url_string)
    val(organism_sci)
  
  output:
    val(fasta_url), emit: reference_fasta_url
    val(gtf_url), emit: reference_gtf_url
    val(gene_annotations_url), emit: gene_annotations_url
    val(reference_source), emit: reference_source
    val(reference_version), emit: reference_version
  
  exec:
    def organisms = [:]
    println "Fetching table from ${annotations_csv_url_string}"
    
    // Check if input is a URL or a local file path
    if (annotations_csv_url_string.startsWith('http://') || annotations_csv_url_string.startsWith('https://')) {
      // For URLs: use toURL() method
      annotations_csv_url_string.toURL().splitEachLine(",") {fields ->
            organisms[fields[1]] = fields
      }
    } else {
      // For local files: use File class
      new File(annotations_csv_url_string).splitEachLine(",") {fields ->
            organisms[fields[1]] = fields
      }
    }
    
    // extract required fields
    organism_key = organism_sci.capitalize().replace("_"," ")
    
    // Check if the organism exists in the table
    if (organisms.containsKey(organism_key)) {
        fasta_url = organisms[organism_key][5]
        gtf_url = organisms[organism_key][6]
        gene_annotations_url = organisms[organism_key][10]
        
        // Convert figshare ndownloader URL to API endpoint
        if (gene_annotations_url != null && gene_annotations_url.contains('figshare.com/ndownloader/files/')) {
            file_id = (gene_annotations_url =~ /.*\/files\/([a-zA-Z0-9]+).*/)[0][1]
            gene_annotations_url = "https://api.figshare.com/v2/file/download/${file_id}"
            println "            Converted figshare URL to API endpoint: ${gene_annotations_url}"
        }
        
        reference_version = organisms[organism_key][3]
        reference_source = organisms[organism_key][4]
        println "Annotation table values parsed for '${organism_key}':"
        println "            Reference Fasta URL: ${fasta_url}"
        println "            Reference GTF URL: ${gtf_url}" 
        println "            Gene Annotations URL: ${gene_annotations_url}"
        println "            Reference Source: ${reference_source}"
        if (reference_source.toLowerCase().contains('ensembl')) {
            println "            Reference Version: ${reference_version}"
        }
    } else {
        println "WARNING: Organism '${organism_key}' not found in annotations table."
        println "Returning null values for all outputs."
        fasta_url = null
        gtf_url = null
        gene_annotations_url = null
        reference_version = null
        reference_source = null
    }
}