process ASSESS_STRANDEDNESS {
  // tag "Dataset-wide"

  input:
    path("infer_out/*") // a collection of infer_experiment stdout files

  output:
    path("result.txt")

  script:
    """
    assess_strandedness.py infer_out ${params.entry_point == 'bam_files' && params.strandedness == 'auto' ? '--strict' : ''}
    """

  stub:
    """
    assess_strandedness.py infer_out ${params.entry_point == 'bam_files' && params.strandedness == 'auto' ? '--strict' : ''}
    echo "unstranded:0.48595" > result.txt # override original results, this is because heavy truncation and genome subsampling can result in an ambiguous strand assignment, which normally is an issue, but should be ignore for stubruns
    """
}